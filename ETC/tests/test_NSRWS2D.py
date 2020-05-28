#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""

from array import array
from random import choice

from hypothesis import given
from hypothesis.strategies import composite, integers, lists

from ETC.NSRWS.x2D import core as cc
from ETC.NSRWS.x2D import etc as cetc
from ETC.NSRWS.x2D.onestep import onestep


@composite
def generate_sequences(draw, elements=[lists, integers]):
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
    seq_x = draw(lists(integers(min_value=1, max_value=1), min_size=3, max_size=10_000))
    seq_y = draw(
        lists(
            integers(min_value=2, max_value=2), min_size=len(seq_x), max_size=len(seq_x)
        )
    )

    return seq_x, seq_y


@given(generate_sequences())
def test_onestep(inputs):
    seq_x, seq_y = inputs
    out_x, out_y, signal = onestep(seq_x, seq_y, order=2, verbose=False, check=False)

    assert len(out_x) < len(seq_x) and len(out_y) < len(seq_y)

    assert max(out_x) > max(seq_x) and max(out_y) > max(seq_y)

    assert min(out_x) >= min(seq_x) and min(out_y) >= min(seq_y)

    assert len(set(out_x) - set(seq_x)) == 1
    assert len(set(out_y) - set(seq_y)) == 1

    assert 1 <= len(set(out_x) ^ set(seq_x)) <= 3
    assert 1 <= len(set(out_y) ^ set(seq_y)) <= 3


@given(generate_sequences())
def test_get_mask_general(inputs):
    seq_x, seq_y = inputs
    seq_x = array("I", seq_x)
    seq_y = array("I", seq_y)
    mask = cc.get_mask_pairs(seq_x, seq_y)

    assert len(mask) == len(seq_x) - 1 == len(seq_y) - 1
    assert set(mask).issubset({0, 1})

    try:
        idx0 = mask.index(0)
        assert seq_x[idx0] == seq_x[idx0 + 1] and seq_y[idx0] == seq_y[idx0 + 1]
    except ValueError:
        pass


@given(generate_sequences_identical())
def test_get_mask_identical(inputs):
    seq_x, seq_y = inputs
    seq_x = array("I", seq_x)
    seq_y = array("I", seq_y)
    mask = cc.get_mask_pairs(seq_x, seq_y)
    mask_rev = cc.get_mask_pairs(seq_x[::-1], seq_y[::-1])

    idx0 = mask.index(0)

    assert len(mask) == len(mask_rev) == len(seq_x) - 1 == len(seq_y) - 1
    assert set(mask).issubset({0, 1})
    assert set(mask_rev).issubset({0, 1})
    assert mask.count(0) == mask_rev.count(0) == (len(seq_x) - 1) // 2
    assert mask.count(1) == mask_rev.count(1) == len(seq_x) // 2

    assert seq_x[idx0] == seq_x[idx0 + 1] and seq_y[idx0] == seq_y[idx0 + 1]


@given(generate_sequences())
def test_substitution(inputs):
    seq_x, seq_y = inputs
    seq_x = array("I", seq_x)
    seq_y = array("I", seq_y)

    sub_value_x = 1 + max(seq_x)
    sub_value_y = 1 + max(seq_y)

    idx = seq_x.index(choice(seq_x[:-1]))
    pair_x = array("I", [seq_x[idx], seq_x[idx + 1]])
    pair_y = array("I", [seq_y[idx], seq_y[idx + 1]])

    out_x, out_y = cc.substitute_pairs(
        seq_x[:], seq_y[:], pair_x, pair_y, sub_value_x, sub_value_y
    )

    assert len(out_x) < len(seq_x) and len(out_y) < len(seq_y)
    assert len(out_x) == len(out_y)
    assert max(out_x) > max(seq_x) and max(out_y) > max(seq_y)


@given(generate_sequences())
def test_truncation(inputs):
    seq_x, seq_y = inputs

    etc_vf = cetc.compute(seq_x, seq_y, order=2, verbose=True, truncate=False)["ETC2D"]
    etc_vt = cetc.compute(seq_x, seq_y, order=2, verbose=True, truncate=True)["ETC2D"]
    etc_cf = cetc.compute(seq_x, seq_y, order=2, verbose=False, truncate=False)["ETC2D"]
    etc_ct = cetc.compute(seq_x, seq_y, order=2, verbose=False, truncate=True)["ETC2D"]

    assert etc_vf == etc_vt == etc_cf == etc_ct


def test_compute_save(tmp_path):

    seq_x = array("I", [2, 4] * 100)
    seq_y = array("I", [1, 3, 5, 7] * 50)
    file = tmp_path / "test.csv"

    etc_vf = cetc.compute_save(seq_x, seq_y, file, order=2, truncate=False)
    assert isinstance(etc_vf, dict)

    etc_vt = cetc.compute_save(seq_x, seq_y, file, order=2, truncate=True)
    assert isinstance(etc_vt, dict)
