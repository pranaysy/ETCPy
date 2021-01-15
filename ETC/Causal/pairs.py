#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""

import ETC
from ETC.NSRWS.x1D import core
from entropy import lziv_complexity as LZ
from gzip import compress

from collections import Counter
from itertools import compress, islice

def check_pair(pair, seq):

    # Create overlapped sliding pairs
    pairs = Counter(zip(*(islice(seq, i, None) for i in range(2)))).keys()

    return pair in pairs


def external_substitution(seq, trajectory):

    # Assign proper type
    seq = ETC.cast(seq)

    # Initialize ETC to 0
    etc = 0

    # Iterate over the given substitution table and substitute
    for pair in trajectory[1:]:  # Skip first entry, not a substitution step

        # Substitute only if the sequence is atleast 2 symbols long
        if len(seq) > 1 and check_pair(tuple(pair.get("window")), seq):

            # Cython function call
            seq = ETC.cast(core.substitute_pairs(seq, pair.get("window"), max(seq) + 1))
            etc += 1

        # If sequence has been fully compressed, stop
        else:
            break

    # Return both etc as well as the sequence, whatever is left of it
    return etc, seq


def ETC_causality(x, y, lengths=False):

    result = {}

    if lengths:
        # Add sequence lengths to results
        result.update({"length_x": len(x), "length_y": len(y)})

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

        etc_y_residual = 0

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

        etc_x_residual = 0

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
    if abs(causal_y_from_x - causal_x_from_y) <= 1:
        result.update({"direction_ETC_penalty": "none_or_mutual"})
    elif causal_y_from_x > causal_x_from_y:
        result.update({"direction_ETC_penalty": "x_causes_y"})
    else:
        result.update({"direction_ETC_penalty": "y_causes_x"})

    # Verbal description of causal direction
    if abs(etc_y_residual - etc_x_residual) <= 1:
        result.update({"direction_ETC_efficacy": "none_or_mutual"})
    elif etc_y_residual < etc_x_residual:
        result.update({"direction_ETC_efficacy": "x_causes_y"})
    else:
        result.update({"direction_ETC_efficacy": "y_causes_x"})

    return result


def LZ_causality(x, y, lengths=False):

    result = {}

    if lengths:
        # Add sequence lengths to results
        result.update({"length_x": len(x), "length_y": len(y)})

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
            "LZ_xy-x": lz_xy - lz_x,
            "LZ_xy-y": lz_xy - lz_y,
            "LZ_xy_norm": lz_xy_norm,
            "LZ_xy-x_norm": lz_xy_norm - lz_x_norm,
            "LZ_xy-y_norm": lz_xy_norm - lz_y_norm,
            "LZ_yx": lz_yx,
            "LZ_yx-y": lz_yx - lz_y,
            "LZ_yx-x": lz_yx - lz_x,
            "LZ_yx_norm": lz_yx_norm,
            "LZ_yx-y_norm": lz_yx_norm - lz_y_norm,
            "LZ_yx-x_norm": lz_yx_norm - lz_x_norm,
        }
    )

    # LZ penalty version
    if abs(result["LZ_xy-y"] - result["LZ_yx-x"]) <= 1:
        LZ_penalty = "none_or_mutual"
    elif result["LZ_xy-y"] > result["LZ_yx-x"]:
        LZ_penalty = "x_causes_y"
    else:
        LZ_penalty = "y_causes_x"

    # LZ penalty version, normalized
    if result["LZ_xy-y_norm"] > result["LZ_yx-x_norm"]:
        LZ_penalty_norm = "x_causes_y"
    elif result["LZ_xy-y_norm"] < result["LZ_yx-x_norm"]:
        LZ_penalty_norm = "y_causes_x"
    else:
        LZ_penalty_norm = "none_or_mutual"

    # LZ efficacy version
    if abs(result["LZ_xy-x"] - result["LZ_yx-y"]) <= 1:
        LZ_efficacy = "none_or_mutual"
    if result["LZ_xy-x"] > result["LZ_yx-y"]:
        LZ_efficacy = "x_causes_y"
    else:
        LZ_efficacy = "y_causes_x"


    # LZ efficacy version, normalized
    if result["LZ_xy-x_norm"] > result["LZ_yx-y_norm"]:
        LZ_efficacy_norm = "x_causes_y"
    elif result["LZ_xy-x_norm"] < result["LZ_yx-y_norm"]:
        LZ_efficacy_norm = "y_causes_x"
    else:
        LZ_efficacy_norm = "none_or_mutual"

    result.update(
        {"direction_LZ_efficacy": LZ_efficacy, "direction_LZ_penalty": LZ_penalty, "direction_LZ_efficacy_norm": LZ_efficacy_norm, "direction_LZ_penalty_norm": LZ_penalty_norm}
    )

    return result

def GZ_causality(x, y):

    gzip = {
        "sz_x" : len(bytearray(x)),
        "sz_x_gz" : len(compress(bytearray(x))),
        "sz_y" : len(bytearray(y)),
        "sz_y_gz" : len(compress(bytearray(y))),
        "sz_xy" : len(bytearray(x)),
        "sz_xy_gz" : len(compress(bytearray(x+y))),
        "sz_yx" : len(bytearray(x)),
        "sz_yx_gz" : len(compress(bytearray(y+x))),
        }

    # LZ efficacy version
    if abs(abs(gzip["sz_yx_gz"] - gzip["sz_y_gz"]) - abs(gzip["sz_xy_gz"] - gzip["sz_x_gz"])) <= 1:
        direction = "none_or_mutual"
    elif (gzip["sz_yx_gz"] - gzip["sz_y_gz"]) > (gzip["sz_xy_gz"] - gzip["sz_x_gz"]):
        direction = "x_causes_y"
    else:
        direction = "y_causes_x"


    gzip.update({"direction_GZ_efficacy":direction})

    return gzip

##-------------------------------------------------------------------------------------#
def LZ_penalty(x, y, lengths=False):

    result = {}

    if lengths:
        # Add sequence lengths to results
        result.update({"length_x": len(x), "length_y": len(y)})

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
            "LZ_xy-y": lz_xy - lz_y,
            "LZ_xy_norm": lz_xy_norm,
            "LZ_xy-y_norm": lz_xy_norm - lz_y_norm,
            "LZ_yx": lz_yx,
            "LZ_yx-x": lz_yx - lz_x,
            "LZ_yx_norm": lz_yx_norm,
            "LZ_yx-x_norm": lz_yx_norm - lz_x_norm,
        }
    )

    # LZ penalty version
    if abs(result["LZ_xy-y"] - result["LZ_yx-x"]) <= 1:
        LZ_penalty = "none_or_mutual"
    elif result["LZ_xy-y"] > result["LZ_yx-x"]:
        LZ_penalty = "x_causes_y"
    else:
        LZ_penalty = "y_causes_x"

    # LZ penalty version, normalized
    if result["LZ_xy-y_norm"] > result["LZ_yx-x_norm"]:
        LZ_penalty_norm = "x_causes_y"
    elif result["LZ_xy-y_norm"] < result["LZ_yx-x_norm"]:
        LZ_penalty_norm = "y_causes_x"
    else:
        LZ_penalty_norm = "none_or_mutual"

    result.update(
        {"direction_LZ_penalty": LZ_penalty, "direction_LZ_penalty_norm": LZ_penalty_norm}
    )

    return result