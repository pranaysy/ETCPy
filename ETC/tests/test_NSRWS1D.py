#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""

from array import array
from random import choice

from hypothesis import given
from hypothesis.strategies import composite, integers, lists

from ETC.NSRWS.x1D import onestep
from ETC.NSRWS.x1D import etc as cetc
from ETC.NSRWS.x1D import core as cc


@composite
def generate_sequence(draw, elements=[lists, integers]):
    """
    Generate a list of integers as sequence input, and an integer >= 2 for order param
    """
    seq = draw(lists(integers(min_value=1, max_value=100), min_size=3, max_size=10_000))

    order = draw(integers(min_value=2, max_value=len(seq) - 1))

    return seq, order


@composite
def generate_sequence_identical(draw, elements=[lists, integers]):
    """
    Generate a list with all identical integers, and an integer >= 2 for order param
    """
    seq = draw(lists(integers(min_value=1, max_value=1), min_size=3, max_size=10_000))

    order = draw(integers(min_value=2, max_value=len(seq) - 1))

    return seq, order


@given(generate_sequence())
def test_onestep(inputs):
    """
    Test the outermost onestep function exposed for direct estimation
    """
    seq, order = inputs
    output1, signal = onestep.onestep(seq, order, verbose=False, check=False)
    output2verbose = onestep.onestep(seq, order, verbose=True, check=False)

    # Substituted sequence should be shorter than input
    assert len(output1) < len(seq)

    # Highest value in substituted sequence should be greater than that in input
    assert max(output1) > max(seq)

    # Smallest value in substituted sequence should be at least as large as that in input
    assert min(output1) >= min(seq)

    # Number of unique symbols in output should be one less than that in input
    assert len(set(output1) - set(seq)) == 1

    # Changing verbosity parameter should not alter the substituted sequence
    assert output1 == output2verbose[0]


@given(generate_sequence())
def test_onestep_invalid(inputs):
    """
    Test the outermost onestep function for invalid input: sequence shorter than order
    """
    seq, order = inputs

    output = onestep.onestep(seq[: order - 1], order, verbose=False, check=False)

    assert output is None


def test_onestep_invalid_str():
    """
    Test the outermost onestep function for invalid input: string input
    """
    output = onestep.onestep("abcdef", 6, verbose=False, check=False)

    assert output is None


@given(generate_sequence_identical())
def test_onestep_identical(inputs):
    """
    Test the outermost onestep function for sequence with identical symbols
    """
    seq, order = inputs

    output = onestep.onestep(seq, order, verbose=False, check=True)

    assert output is None


@given(generate_sequence())
def test_onestep_pairs_vs_windows(inputs):
    """
    Test the parity of output for NSRWS with order=2 & NSRPS for random sequences
    """
    seq, _ = inputs
    seq = array("I", seq)

    # Pairs (NSRPS)
    out_pairs = onestep._onestep_pairs(seq[:], verbose=True)

    # Order = 2 (NSRWS)
    out_windows = onestep._onestep_windows(seq[:], 2, verbose=True)

    # Check equality of all but last (timings) part of output
    for n in range(4):
        assert out_pairs[n] == out_windows[n]


@given(generate_sequence_identical())
def test_onestep_pairs_vs_windows_identical(inputs):
    """
    Test the parity of output for NSRWS with order=2 & NSRPS for all identical sequence
    """
    seq, _ = inputs
    seq = array("I", seq)

    # Pairs (NSRPS)
    out_pairs = onestep._onestep_pairs(seq[:], verbose=True)

    # Order = 2 (NSRWS)
    out_windows = onestep._onestep_windows(seq[:], 2, verbose=True)

    # Check equality of all but last (timings) part of output
    for n in range(4):
        assert out_pairs[n] == out_windows[n]


@given(generate_sequence())
def test_get_mask_general(inputs):
    """
    Test the get_mask function for random sequences and orders
    """
    seq, order = inputs
    seq = array("I", seq)

    # Get mask depending on order
    if order == 2:
        mask = cc.get_mask_pairs(seq)
    else:
        mask = cc.get_mask_windows(seq, order)

    # Mask should be precisely shorter than input sequence
    assert len(mask) == len(seq) - order + 1

    # Mask should only contain 0s and 1s
    assert set(mask).issubset({0, 1})

    # First element must be 1
    assert mask[0] == 1

    # If mask contains a 0, then that position in the sequence indicates an overlap
    try:
        idx0 = mask.index(0)
        # Check if consecutive elements equal where 0 found in mask
        assert seq[idx0] == seq[idx0 + order - 1]
    except ValueError:
        pass


@given(generate_sequence_identical())
def test_get_mask_identical(inputs):
    """
    Test the get_mask function for sequences with identical symbols
    """
    seq, order = inputs
    seq = array("I", seq)

    # Get mask depending on order from left-to-right and reversed sequence
    if order == 2:
        mask = cc.get_mask_pairs(seq)
        mask_rev = cc.get_mask_pairs(seq[::-1])
    else:
        mask = cc.get_mask_windows(seq, order)
        mask_rev = cc.get_mask_windows(seq[::-1], order)

    # Find zeroes for they must be present
    idx0 = mask.index(0)

    # Both masks should be precisely shorter than input sequence
    assert len(mask) == len(mask_rev) == len(seq) - order + 1

    # Both masks should only contain 0s and 1s
    assert set(mask).issubset({0, 1})
    assert set(mask_rev).issubset({0, 1})

    # First element has to be a 1
    assert mask[0] == 1

    # Check if consecutive elements equal where 0 found in mask
    assert seq[idx0] == seq[idx0 + order - 1]


def test_mask_and_count():
    """
    Test the function for applying mask and counting frequent windows
    """
    seq = (1, 2, 3, 4, 5, 6, 7)
    mask = (1, 0, 0, 1, 1)
    assert onestep._mask_and_count(seq, mask, 3) == (array("I", (1, 2, 3)), 1)

    seq = (1, 2, 3, 4, 5, 6, 7)
    mask = (1, 1, 1, 1, 1)
    assert onestep._mask_and_count(seq, mask, 3) == (array("I", (1, 2, 3)), 1)

    seq = (1, 2, 3, 4, 5, 6, 7)
    mask = (0, 1, 1, 1, 1)
    assert onestep._mask_and_count(seq, mask, 3) == (array("I", (2, 3, 4)), 1)

    seq = (1, 1, 1, 1, 1, 2, 1)
    mask = (1, 0, 0, 1, 1)
    assert onestep._mask_and_count(seq, mask, 3) == (array("I", (1, 1, 1)), 1)


@given(generate_sequence())
def test_substitution(inputs):
    """
    Test the substitution step for random sequences
    """
    seq, order = inputs
    seq = array("I", seq)

    # Get value to substitute
    sub_value = 1 + max(seq)

    # Pick a random pair for substitution
    idx = seq.index(choice(seq[:-1]))
    pair = array("I", [seq[idx], seq[idx + 1]])

    # Substitute the pair using both functions
    out1 = cc.substitute_pairs(seq[:], pair, sub_value)
    out2 = cc.substitute_windows(seq[:], 2, pair, sub_value)

    # The 2 outputs should be equal
    assert out1 == out2

    # The length of the substituted sequence should be less than the input sequence
    assert len(out1) < len(seq)

    # The highest value in the substituted sequence should be more than that in the input sequence
    assert max(out1) > max(seq)

    # The highest value in the substitute sequence should match the provided value
    assert max(out1) == sub_value


@given(generate_sequence())
def test_truncation(inputs):
    """
    Test ETC estimation from all 4 methods based on verbosity and truncation
    """
    seq, order = inputs

    etc_vf = cetc.compute(seq, order, verbose=True, truncate=False)["ETC1D"]
    etc_vt = cetc.compute(seq, order, verbose=True, truncate=True)["ETC1D"]
    etc_cf = cetc.compute(seq, order, verbose=False, truncate=False)["ETC1D"]
    etc_ct = cetc.compute(seq, order, verbose=False, truncate=True)["ETC1D"]

    # All 4 estimates should be identical
    assert etc_vf == etc_vt == etc_cf == etc_ct


def test_compute_save(tmp_path):
    """
    Test ETC estimation with write-to-disk functionality
    """
    seq = array("I", [2, 4] * 100)

    # Temporary file (Path object) for use
    file = tmp_path / "test.csv"

    # Test without truncation
    etc_vf = cetc.compute_save(seq, file, order=2, truncate=False)
    assert isinstance(etc_vf, dict)

    # Test with truncation
    etc_vt = cetc.compute_save(seq, file, order=2, truncate=True)
    assert isinstance(etc_vt, dict)

    # Values should be same of course
    assert etc_vf["ETC1D"] == etc_vt["ETC1D"]
