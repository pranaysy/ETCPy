#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""

from array import array
from random import choice

from hypothesis import given
from hypothesis.strategies import composite, integers, lists

from ETC.NSRWS.x2D import onestep
from ETC.NSRWS.x2D import core as cc
from ETC.NSRWS.x2D import etc as cetc


@composite
def generate_sequences(draw, elements=[lists, integers]):
    """
    Generate 2 lists of integers as sequence input with equal lengths
    """
    seq_x = draw(
        lists(integers(min_value=1, max_value=100), min_size=3, max_size=10_000)
    )
    seq_y = draw(
        lists(
            integers(min_value=1, max_value=100),
            min_size=len(seq_x),
            max_size=len(seq_x),
        )
    )

    return seq_x, seq_y


@composite
def generate_sequences_identical(draw, elements=[lists, integers]):
    """
    Generate 2 lists with all identical integers with equal lengths
    """
    seq_x = draw(lists(integers(min_value=1, max_value=1), min_size=3, max_size=10_000))
    seq_y = draw(
        lists(
            integers(min_value=2, max_value=2), min_size=len(seq_x), max_size=len(seq_x)
        )
    )

    return seq_x, seq_y


@given(generate_sequences())
def test_onestep(inputs):
    """
    Test the outermost onestep function exposed for direct estimation
    """
    seq_x, seq_y = inputs

    out_x1, out_y1, _ = onestep.onestep(
        seq_x, seq_y, order=2, verbose=False, check=False
    )
    outputs = onestep.onestep(seq_x, seq_y, order=2, verbose=True, check=False)
    out_x2, out_y2 = outputs[0], outputs[1]

    # Substituted sequence should be shorter than input
    assert len(out_x1) < len(seq_x) and len(out_y1) < len(seq_y)

    # Highest value in substituted sequence should be greater than that in input
    assert max(out_x1) > max(seq_x) and max(out_y1) > max(seq_y)

    # Smallest value in substituted sequence should be at least as large as that in input
    assert min(out_x1) >= min(seq_x) and min(out_y1) >= min(seq_y)

    # Number of unique symbols in output should be one less than that in input
    assert len(set(out_x1) - set(seq_x)) == 1
    assert len(set(out_y1) - set(seq_y)) == 1

    # Number of symbols in output that are not in input should be between 1 and 3
    assert 1 <= len(set(out_x1) ^ set(seq_x)) <= 3
    assert 1 <= len(set(out_y1) ^ set(seq_y)) <= 3

    # Changing verbosity parameter should not alter the substituted sequence
    assert out_x1 == out_x2 and out_y1 == out_y2


@given(generate_sequences())
def test_onestep_unequal(inputs):
    """
    Test the outermost onestep function for invalid input: sequences shorter than order
    """
    seq_x, seq_y = inputs

    output = onestep.onestep(seq_x, seq_y[:-11], order=2, verbose=False, check=False)
    assert output is None


@given(generate_sequences())
def test_onestep_invalid(inputs):
    """
    Test the outermost onestep function for invalid input: sequences shorter than order
    """
    seq_x, seq_y = inputs

    output = onestep.onestep(seq_x[:1], seq_y[:1], order=2, verbose=False, check=False)
    assert output is None


def test_onestep_invalid_str():
    """
    Test the outermost onestep function for invalid input: string inputs
    """
    output = onestep.onestep(
        [1, 2, 3, 4, 5, 6], "abcdef", 2, verbose=False, check=False
    )
    assert output is None

    output = onestep.onestep(
        "abcdef", [1, 2, 3, 4, 5, 6], 2, verbose=False, check=False
    )
    assert output is None

    output = onestep.onestep("abcdef", "abcdef", 2, verbose=False, check=False)
    assert output is None


@given(generate_sequences_identical())
def test_onestep_identical(inputs):
    """
    Test the outermost onestep function for sequence with identical symbols
    """
    seq_x, seq_y = inputs

    output = onestep.onestep(seq_x, seq_y, order=2, verbose=False, check=True)

    assert output is None


@given(generate_sequences())
def test_get_mask_general(inputs):
    """
    Test the get_mask function for random sequences and orders
    """
    seq_x, seq_y = inputs
    seq_x = array("I", seq_x)
    seq_y = array("I", seq_y)

    # Get mask
    mask = cc.get_mask_pairs(seq_x, seq_y)

    # Mask should be precisely shorter than input sequence
    assert len(mask) == len(seq_x) - 1 == len(seq_y) - 1

    # Mask should only contain 0s and 1s
    assert set(mask).issubset({0, 1})

    # First element must be 1
    assert mask[0] == 1

    # If mask contains a 0, then that position in the sequence indicates an overlap
    try:
        idx0 = mask.index(0)
        # Check if consecutive elements equal where 0 found in mask
        assert seq_x[idx0] == seq_x[idx0 + 1] and seq_y[idx0] == seq_y[idx0 + 1]
    except ValueError:
        pass


@given(generate_sequences_identical())
def test_get_mask_identical(inputs):
    """
    Test the get_mask function for sequences with identical symbols
    """
    seq_x, seq_y = inputs
    seq_x = array("I", seq_x)
    seq_y = array("I", seq_y)

    # Get mask from left-to-right and reversed sequence
    mask = cc.get_mask_pairs(seq_x, seq_y)
    mask_rev = cc.get_mask_pairs(seq_x[::-1], seq_y[::-1])

    # Find zeroes for they must be present
    idx0 = mask.index(0)

    # Both masks should be precisely shorter than input sequence
    assert len(mask) == len(mask_rev) == len(seq_x) - 1 == len(seq_y) - 1

    # Both masks should only contain 0s and 1s
    assert set(mask).issubset({0, 1})
    assert set(mask_rev).issubset({0, 1})

    # Both masks should have an exact number of zeros corresponding to overlaps
    assert mask.count(0) == mask_rev.count(0) == (len(seq_x) - 1) // 2
    assert mask.count(1) == mask_rev.count(1) == len(seq_x) // 2

    # Check if consecutive elements equal where 0 found in mask
    assert seq_x[idx0] == seq_x[idx0 + 1] and seq_y[idx0] == seq_y[idx0 + 1]


def test_mask_and_count():
    """
    Test the function for applying mask and counting frequent windows
    """
    seq_x = (1, 2, 3, 4, 5, 6, 7)
    seq_y = (3, 4, 5, 6, 7, 8, 9)
    mask = (1, 0, 0, 1, 1)

    freq_pair_x, freq_pair_y, count = onestep._mask_and_count(seq_x, seq_y, mask, 2)
    assert (
        freq_pair_x == array("I", (1, 2))
        and freq_pair_y == array("I", (3, 4))
        and count == 1
    )

    mask = (1, 1, 1, 1, 1)

    freq_pair_x, freq_pair_y, count = onestep._mask_and_count(seq_x, seq_y, mask, 2)
    assert (
        freq_pair_x == array("I", (1, 2))
        and freq_pair_y == array("I", (3, 4))
        and count == 1
    )

    mask = (0, 1, 1, 1, 1)
    freq_pair_x, freq_pair_y, count = onestep._mask_and_count(seq_x, seq_y, mask, 2)
    assert (
        freq_pair_x == array("I", (2, 3))
        and freq_pair_y == array("I", (4, 5))
        and count == 1
    )


@given(generate_sequences())
def test_substitution(inputs):
    """
    Test the substitution step for random sequences
    """
    seq_x, seq_y = inputs
    seq_x = array("I", seq_x)
    seq_y = array("I", seq_y)

    # Get values to substitute
    sub_value_x = 1 + max(seq_x)
    sub_value_y = 1 + max(seq_y)

    # Pick a random pair for substitution
    idx = seq_x.index(choice(seq_x[:-1]))
    pair_x = array("I", [seq_x[idx], seq_x[idx + 1]])
    pair_y = array("I", [seq_y[idx], seq_y[idx + 1]])

    # Substitute the pairs
    out_x, out_y = cc.substitute_pairs(
        seq_x[:], seq_y[:], pair_x, pair_y, sub_value_x, sub_value_y
    )

    # The length of the substituted sequence should be less than the input sequence
    assert len(out_x) < len(seq_x) and len(out_y) < len(seq_y)

    # The lengths of the 2 substituted sequences should be identical
    assert len(out_x) == len(out_y)

    # The highest value in the substituted sequence should be more than that in the input sequence
    assert max(out_x) > max(seq_x) and max(out_y) > max(seq_y)

    # The highest value in the substitute sequence should match the provided value
    assert max(out_x) == sub_value_x and max(out_y) == sub_value_y


@given(generate_sequences())
def test_truncation(inputs):
    """
    Test ETC estimation from all 4 methods based on verbosity and truncation
    """
    seq_x, seq_y = inputs

    etc_vf = cetc.compute(seq_x, seq_y, order=2, verbose=True, truncate=False)["ETC2D"]
    etc_vt = cetc.compute(seq_x, seq_y, order=2, verbose=True, truncate=True)["ETC2D"]
    etc_cf = cetc.compute(seq_x, seq_y, order=2, verbose=False, truncate=False)["ETC2D"]
    etc_ct = cetc.compute(seq_x, seq_y, order=2, verbose=False, truncate=True)["ETC2D"]

    # All 4 estimates should be identical
    assert etc_vf == etc_vt == etc_cf == etc_ct


def test_compute_save(tmp_path):
    """
    Test ETC estimation with write-to-disk functionality
    """
    seq_x = array("I", [2, 4] * 100)
    seq_y = array("I", [1, 3, 5, 7] * 50)

    # Temporary file (Path object) for use
    file = tmp_path / "test.csv"

    # Test without truncation
    etc_vf = cetc.compute_save(seq_x, seq_y, file, order=2, truncate=False)
    assert isinstance(etc_vf, dict)

    # Test with truncation
    etc_vt = cetc.compute_save(seq_x, seq_y, file, order=2, truncate=True)
    assert isinstance(etc_vt, dict)

    # Values should be same of course
    assert etc_vf["ETC2D"] == etc_vt["ETC2D"]
