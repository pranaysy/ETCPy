#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module contains functions for causal discovery from a pair of symbolic sequences
using compression-complexity measures such as ETC and LZ complexity.

Three functions are exposed for application on array.array data types:
    ETC_causality: Computes estimates from 2 models - penalty & efficacy using ETC
    LZ_causality: Computes estimates from 1 model - penalty using LZ
    CCM_causality: Wrapper for both the above

Each model returns a causal direction and strength of evidence in favor of that direction

@author: Pranay S. Yadav
"""
# Import libraries
import ETC
from ETC.NSRWS.x1D import core
from ETC.LZ76.lzc import compute_complexity as LZ

# from entropy import lziv_complexity as LZ

# from collections import Counter
from itertools import islice
from hashlib import blake2b

# Check if in input pair is present in a sequence
def _check_pair(pair, seq):
    """
    Check whether a pair is present in an input sequence

    Uses itertools tricks to create staggered / sliding sequence tuples and gets counts
    using Counter from collections. The keys to

    Parameters
    ----------
    pair : tuple
        Pair of elements to check.
    seq : array.array, tuple
        Discrete symbolic sequence containing possible pair.

    Returns
    -------
    bool
        True if pair present in seq.

    """

    ##### OLD METHOD using Counter #####################################################
    # Create overlapped sliding pairs, count and get keys
    # pairs = Counter(zip(*(islice(seq, i, None) for i in range(2)))).keys()

    ##### NEW METHOD using set #########################################################
    # Convert seq to tuple, create overlapped sliding pairs, get unique tuples with set
    pairs = set(zip(*(islice(tuple(seq), i, None) for i in range(2))))

    # Return a boolean
    return pair in pairs


def _external_substitution(seq, trajectory):
    """
    Carry out pair substitution for a sequence given a trajectory / successive sequence
    of pairs obtained externally (eg from compression of another sequence)

    For each pair present in the trajectory (except the first one), check if it is present
    in the sequence, and if present substitute it as per standard ETC routines. Conditions
    for stopping, any one is sufficient:
        1. Sequence gets fully compressed (substituted) so that it's length is equal to 1
        2. Pairs get fully exhausted, no more substitution possible
        3. Pairs not present in the sequence

    Uses _check_pair for testing existence of pair in seq & Cython calls for substitution

    Parameters
    ----------
    seq : array.array
        Discrete symbolic sequence containing 32-bit unsigned integers.
    trajectory : list
        List of dictionaries corresponding to each step of NSRWS run during
        estimation of ETC for the given sequence.

    Returns
    -------
    etc : int
        ETC for cross-compression of seq using trajectory.
    seq : array.array
        Residual discrete symbolic sequence left over after compression using trajectory.

    """
    # Assign proper type
    seq = ETC.cast(seq)

    # Initialize ETC to 0
    etc = 0

    # Iterate over the given substitution table and substitute
    for step in trajectory[1:]:  # Skip first entry, not a substitution step

        pair = step.get("window")

        # Substitute only if the sequence is atleast 2 symbols long
        if len(seq) > 1 and _check_pair(tuple(pair), seq):

            # Cython function call
            seq = ETC.cast(core.substitute_pairs(seq, pair, max(seq) + 1))
            etc += 1

        # If sequence has been fully compressed, stop
        else:
            break

    # Return both etc as well as the sequence, whatever is left of it
    return etc, seq


def _ETC_residual(residual_sequence):
    """
    Compute ETC of the residual sequence

    Essentially just the regular ETC, but with additional checks

    Parameters
    ----------
    residual_sequence : array.array
        Residual discrete symbolic sequence left over after compression using trajectory.

    Returns
    -------
    etc : int
        ETC of the input residual_sequence.

    """
    # If the y's residual sequence is long enough, then compress it, get causal estimate
    if len(residual_sequence) > 1:

        # Compress
        return ETC.compute_1D(residual_sequence).get("ETC1D")

    # Already compressed, no residual left
    return 0


def ETC_causality(x, y, penalty_threshold=1, efficacy_tolerance=0, lengths=True):
    """
    Causal discovery and estimation using ETC for a pair of discrete symbolic sequences

    Inference of causal discovery as well as estimation of causal strength of interaction
    from a pair of sequences based on two different, although related models - penalty
    and efficacy.

    Parameters
    ----------
    x, y : array.array, list or tuple
        Sequence of integers.
    penalty_threshold : int, optional, non-negative
        Threshold for difference in causal estimates based on penalty. The default is 1.
    efficacy_tolerance : float, optional, non-negative, between 0 and 1
        Tolerance for difference in causal estimates based on efficacy. The default is 0.
    lengths : bool, optional
        Whether to add sequence lengths to output dict or not. The default is True.
        Useful when both LZ_causality and ETC_causality are called, avoids duplication

    Returns
    -------
    result : dict
        Various estimates from both penalty (ETCP) and efficacy models (ETCE). Includes
        1D-ETC estimates for both sequences; direction & strengths of causal interaction

    """
    # If either was not successfully converted to array, break
    if not (ETC.check.arraytype(x) and ETC.check.arraytype(y)):

        # Convert inputs to arrays
        x = ETC.cast(x)
        y = ETC.cast(y)

        # If unsuccessful, break
        if (x is None) or (y is None):
            return None

    # Initialize output dictionary
    result = {}

    # Add sequence lengths if requested
    if lengths:
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
    etc_y_given_x, y_residual = _external_substitution(y, out_x.get("Trajectory"))
    etc_x_given_y, x_residual = _external_substitution(x, out_y.get("Trajectory"))

    # Add to results
    result.update(
        {"ETC_x_given_Gy": etc_x_given_y, "ETC_y_given_Gx": etc_y_given_x,}
    )

    # Compute ETC for the 2 residuals
    etc_y_residual = _ETC_residual(y_residual)
    etc_x_residual = _ETC_residual(x_residual)

    # Compute penalty estimates in each direction
    penalty_x_to_y = etc_y_given_x + etc_y_residual - etc_y
    penalty_y_to_x = etc_x_given_y + etc_x_residual - etc_x

    # Compute efficacy estimates in each direction, with default normalization
    # Only for non-empty / non-zero residuals
    if etc_y_residual != 0:
        ETC_y_residual_norm = etc_y_residual / (len(y_residual) - 1)
        efficacy_x_to_y = etc_y_residual / (len(y_residual) - 1)
    else:
        ETC_y_residual_norm = 0
        efficacy_x_to_y = 0
    if etc_x_residual != 0:
        ETC_x_residual_norm = etc_x_residual / (len(x_residual) - 1)
        efficacy_y_to_x = etc_x_residual / (len(x_residual) - 1)
    else:
        ETC_x_residual_norm = 0
        efficacy_y_to_x = 0

    # Add to results: residuals and penalty estimates
    result.update(
        {
            "length_x_residual": len(x_residual),
            "length_y_residual": len(y_residual),
            "ETC_x_residual": etc_x_residual,
            "ETC_y_residual": etc_y_residual,
            "ETC_x_residual_norm": ETC_x_residual_norm,
            "ETC_y_residual_norm": ETC_y_residual_norm,
            "ETCP_x_to_y": penalty_x_to_y,
            "ETCP_y_to_x": penalty_y_to_x,
            "ETCP_strength": abs(penalty_y_to_x - penalty_x_to_y),
            "ETCP_threshold": penalty_threshold,
        }
    )

    # Textual description of causal direction based on penalty
    if abs(penalty_y_to_x - penalty_x_to_y) <= penalty_threshold:
        result.update({"ETCP_direction": "n_or_m", "ETCP_cause": "n_or_m"})

    elif penalty_y_to_x > penalty_x_to_y:
        result.update({"ETCP_direction": "x_causes_y", "ETCP_cause": "x"})

    else:
        result.update({"ETCP_direction": "y_causes_x", "ETCP_cause": "y"})

    # Add to results: efficacy estimates
    result.update(
        {
            "ETCE_x_to_y": efficacy_x_to_y,
            "ETCE_y_to_x": efficacy_y_to_x,
            "ETCE_strength": abs(efficacy_x_to_y - efficacy_y_to_x),
            "ETCE_tolerance": efficacy_tolerance,
        }
    )

    # Verbal description of causal direction based on efficacy
    if abs(efficacy_x_to_y - efficacy_y_to_x) <= efficacy_tolerance:
        result.update({"ETCE_direction": "n_or_m", "ETCE_cause": "n_or_m"})

    elif efficacy_x_to_y > efficacy_y_to_x:
        result.update({"ETCE_direction": "x_causes_y", "ETCE_cause": "x"})

    else:
        result.update({"ETCE_direction": "y_causes_x", "ETCE_cause": "y"})

    return result


# def _convert_to_list(seq):
#     """
#     Convert a sequence to a list

#     LZ function from entropy package needs a list or a string as input

#     Parameters
#     ----------
#     seq : array.array, tuple, list, str
#         Discrete symbolic sequence for compression using LZ.

#     Returns
#     -------
#     list
#         Input sequence cast into list.

#     """
#     if not isinstance(seq, list):
#         try:
#             return list(seq)
#         except:  # catch-all, bad practice, I know
#             print("> ERROR: Could not convert input sequence to list!")
#             return seq
#     return seq


def LZ_causality(x, y, penalty_threshold=1, lengths=True):
    """
    Causal discovery and estimation using LZ for a pair of discrete symbolic sequences

    Inference of causal discovery as well as estimation of causal strength of interaction
    from a pair of sequences based on a penalty approach similar to conditional entropy.

    Original idea developed by Nithin Nagaraj based on preliminary experiments

    Parameters
    ----------
    x, y : array.array, list or tuple
        Sequence of integers.
    penalty_threshold : int, optional, non-negative
        Threshold for difference in causal estimates based on penalty. The default is 1.
    lengths : bool, optional
        Whether to add sequence lengths to output dict or not. The default is True.
        Useful when both LZ_causality and ETC_causality are called, avoids duplication

    Returns
    -------
    result : dict
        Various estimates from the penalty (LZP) model. Includes LZ complexity estimates
        for both sequences; direction & strengths of causal interaction

    """
    # If either was not successfully converted to array, break
    if not (ETC.check.arraytype(x) and ETC.check.arraytype(y)):

        # Convert inputs to arrays
        x = ETC.cast(x)
        y = ETC.cast(y)

        # If unsuccessful, break
        if (x is None) or (y is None):
            return None

    result = {}

    # Add sequence lengths if requested
    if lengths:
        result.update({"length_x": len(x), "length_y": len(y)})

    # Compute LZ for the 2 sequences
    LZ_x = LZ(x)
    LZ_y = LZ(y)

    # Compute LZ after concatenation in both directions
    LZ_xy = LZ(x + y)
    LZ_yx = LZ(y + x)
    LZ_concat_mean = 0.5 * (LZ_xy + LZ_yx)
    LZ_concat_diff = abs(LZ_xy - LZ_yx)

    # Conditioning
    LZ_x_given_y = LZ_xy - LZ_y
    LZ_y_given_x = LZ_yx - LZ_x

    result.update(
        {
            "LZ_x": LZ_x,
            "LZ_y": LZ_y,
            "LZ_xy": LZ_xy,
            "LZ_yx": LZ_yx,
            "LZ_concat_mean": LZ_concat_mean,
            "LZ_concat_diff": LZ_concat_diff,
            "LZP_x_to_y": LZ_y_given_x,  # x->y is y|x (subtract x)
            "LZP_y_to_x": LZ_x_given_y,  # y->x is x|y (subtract y)
            "LZP_strength": abs(LZ_x_given_y - LZ_y_given_x),
            "LZP_threshold": penalty_threshold,
        }
    )

    # LZ penalty version
    if abs(LZ_x_given_y - LZ_y_given_x) <= penalty_threshold:
        result.update({"LZP_direction": "n_or_m", "LZP_cause": "n_or_m"})

    elif LZ_x_given_y > LZ_y_given_x:
        result.update({"LZP_direction": "x_causes_y", "LZP_cause": "x"})

    else:
        result.update({"LZP_direction": "y_causes_x", "LZP_cause": "y"})

    return result


def CCM_causality(x, y, penalty_threshold=1, efficacy_tolerance=0, hashes=False):
    """
    Wrapper around ETC_causality and LZ_causality, calls both and returns merged estimates
    for causal discovery based on ETCP, ETCE and LZP models

    Parameters
    ----------
    x, y : array.array, list or tuple
        Sequence of integers.
    penalty_threshold : int, optional, non-negative
        Threshold for difference in causal estimates based on penalty. The default is 1.
    efficacy_tolerance : float, optional, non-negative, between 0 and 1
        Tolerance for difference in causal estimates based on efficacy. The default is 0.

    Returns
    -------
    result : dict
        Various estimates from both penalty (ETCP as well as LZP) and efficacy models
        (ETCE). Includes 1D-ETC as well as LZ complexity estimates for both sequences;
        direction & strengths of causal interaction

    """
    result = {}

    # Add hashes of sequences
    if hashes:
        result.update(
            {
                "x_hash": blake2b(bytearray(x), digest_size=32).hexdigest(),
                "y_hash": blake2b(bytearray(y), digest_size=32).hexdigest(),
                "hash_algorithm": "blake2b",
                "digest_size": 32,
            }
        )

    result.update(
        ETC_causality(
            x,
            y,
            penalty_threshold=penalty_threshold,
            efficacy_tolerance=efficacy_tolerance,
            lengths=True,
        )
    )
    result.update(
        LZ_causality(x, y, penalty_threshold=penalty_threshold, lengths=False,)
    )

    if result["ETCP_cause"] == result["ETCE_cause"] == result["LZP_cause"]:
        result.update({"Consensus": True})
    else:
        result.update({"Consensus": False})

    return result
