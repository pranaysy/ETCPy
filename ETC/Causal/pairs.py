#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""

import ETC
from ETC.NSRWS.x1D import core
from entropy import lziv_complexity as LZ


def external_substitution(seq, trajectory):

    # Assign proper type
    seq = ETC.cast(seq)

    # Initialize ETC to 0
    etc = 0

    # Iterate over the given substitution table and substitute
    for pair in trajectory[1:]:  # Skip first entry, not a substitution step

        # Substitute only if the sequence is atleast 2 symbols long
        if len(seq) > 1:

            # Cython function call
            seq = ETC.cast(core.substitute_pairs(seq, pair.get("window"), max(seq) + 1))
            etc += 1

        # If sequence has been fully compressed, stop
        else:
            break

    # Return both etc as well as the sequence, whatever is left of it
    return etc, seq


def ETC_causality(x, y):

    # Add sequence lengths to results
    result = {"length_x": len(x), "length_y": len(y)}

    # Compute ETC for the 2 sequences
    out_x = ETC.compute_1D(x, order=2, verbose=True, truncate=False)
    out_y = ETC.compute_1D(y, order=2, verbose=True, truncate=False)

    # Store the estimates separately
    etc_x = out_x.get("ETC1D")
    etc_y = out_y.get("ETC1D")

    # Add them to the results
    result.update({"ETC_x": etc_x, "ETC_y": etc_y})

    # Use each other's substitution tables for compressing the other
    etc_y_given_x, y_residual = external_substitution(y, out_x.get("Trajectory"))
    etc_x_given_y, x_residual = external_substitution(x, out_y.get("Trajectory"))

    # Add to results
    result.update(
        {
            "ETC_x_compressedby_ytable": etc_x_given_y,
            "ETC_y_compressedby_xtable": etc_y_given_x,
        }
    )

    # If the y's residual sequence is long enough, then compress it, get causal estimate
    if len(y_residual) > 1:

        # Compress
        etc_y_residual = ETC.compute_1D(y_residual).get("ETC1D")

        # Add this residual's ETC to get total and compute formulated causal estimate
        causal_y_from_x = etc_y_given_x + etc_y_residual - etc_y

        # Add to results
        result.update(
            {
                "length_y_residual": len(y_residual),
                "ETC_y_residual": etc_y_residual,
                "y_caused_by_x": causal_y_from_x,
            }
        )
    # Else, don't compress and get causal estimate
    else:
        # As formulated
        causal_y_from_x = etc_y_given_x - etc_y

        # Add to results
        result.update(
            {
                "length_y_residual": len(y_residual),
                "ETC_y_residual": 0,
                "y_caused_by_x": causal_y_from_x,
            }
        )

    # If the x's residual sequence is long enough, then compress it, get causal estimate
    if len(x_residual) > 1:

        # Compress
        etc_x_residual = ETC.compute_1D(x_residual).get("ETC1D")

        # Add this residual's ETC to get total and compute formulated causal estimate
        causal_x_from_y = etc_x_given_y + etc_x_residual - etc_x

        # Add to results
        result.update(
            {
                "length_x_residual": len(x_residual),
                "ETC_x_residual": etc_x_residual,
                "x_caused_by_y": causal_x_from_y,
            }
        )

    # Else, don't compress and get causal estimate
    else:

        # As formulated
        causal_x_from_y = etc_x_given_y - etc_x

        # Add to results
        result.update(
            {
                "length_x_residual": len(x_residual),
                "ETC_x_residual": 0,
                "x_caused_by_y": causal_x_from_y,
            }
        )

    # Verbal description of causal direction
    if causal_y_from_x == causal_x_from_y:
        result.update({"direction": "none_or_mutual"})
    elif causal_y_from_x > causal_x_from_y:
        result.update({"direction": "x_causes_y"})
    else:
        result.update({"direction": "y_causes_x"})

    return result


def LZ_causality(x, y):

    # Add sequence lengths to results
    result = {"length_x": len(x), "length_y": len(y)}

    # Compute ETC for the 2 sequences
    lz_x = LZ(x, normalize=False)
    lz_x_norm = LZ(x, normalize=True)

    lz_y = LZ(y, normalize=False)
    lz_y_norm = LZ(y, normalize=True)

    result.update(
        {"LZ_x": lz_x, "LZ_x_norm": lz_x_norm, "LZ_y": lz_y, "LZ_y_norm": lz_y_norm}
    )

    lz_xy = LZ(x + y, normalize=False)
    lz_xy_norm = LZ(x + y, normalize=True)

    lz_yx = LZ(y + x, normalize=False)
    lz_yx_norm = LZ(y + x, normalize=True)

    result.update(
        {
            "LZ_xy": lz_xy,
            "LZ_xy_norm": lz_xy_norm,
            "LZ_yx": lz_yx,
            "LZ_yx_norm": lz_yx_norm,
        }
    )

    return result
