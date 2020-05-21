#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""

from array import array
from random import choice

from hypothesis import given
from hypothesis.strategies import composite, integers, lists

from ETC.NSRWS.x1D.compute_one_step import one_step
from ETC.NSRWS.x1D import compute_etc as cetc
from ETC.NSRWS.x1D import compute_core as cc


@composite
def generate_sequence(draw, elements=[lists, integers]):
    seq = draw(lists(integers(min_value=1, max_value=100), min_size=3, max_size=10_000))

    order = draw(integers(min_value=2, max_value=len(seq)-1))

    return seq, order

@composite
def generate_sequence_identical(draw, elements=[lists, integers]):
    seq = draw(lists(integers(min_value=1, max_value=1), min_size=3, max_size=10_000))
    order = draw(integers(min_value=2, max_value=len(seq)-1))

    return seq, order


@given(generate_sequence())
def test_one_step(inputs):
    x, order = inputs
    output1, signal = one_step(x, order, verbose=False, check=False)
    output2 = one_step(x, order, verbose=True, check=False)

    assert len(output1) < len(x)
    assert max(output1) > max(x)
    assert min(output1) >= min(x)
    assert len(set(output1) - set(x)) == 1

    assert output1 == output2[0]


@given(generate_sequence_identical())
def test_one_step_equality(inputs):

    seq, order = inputs
    output = one_step(seq, order, verbose=False, check=True)

    assert output is None

@given(generate_sequence())
def test_get_mask_general(inputs):
    x, order = inputs
    x = array("I", x)

    if order == 2:
        mask = cc.get_mask_pairs(x)[:-1]
    else:
        mask = cc.get_mask_windows(x, order)[:-(order-1)]


    assert len(mask) < len(x)
    assert set(mask).issubset({0, 1})

    try:
        idx0 = mask.index(0)
        assert x[idx0] == x[idx0 + 1]
    except ValueError:
        pass

@given(generate_sequence_identical())
def test_get_mask_identical(inputs):
    x, order = inputs
    x = array("I", x)

    if order == 2:
        mask = cc.get_mask_pairs(x)[:-1]
        mask_rev = cc.get_mask_pairs(x[::-1])[:-1]
    else:
        mask = cc.get_mask_windows(x, order)[:-(order-1)]
        mask_rev = cc.get_mask_windows(x[::-1], order)[:-(order-1)]

    idx0 = mask.index(0)

    assert len(mask) == len(mask_rev) == len(x) - order + 1
    assert set(mask).issubset({0, 1})
    assert set(mask_rev).issubset({0, 1})

    assert x[idx0] == x[idx0 + 1]


@given(generate_sequence())
def test_substitution(inputs):
    x, order = inputs
    x = array("I", x)

    sub_value = 1 + max(x)

    idx = x.index(choice(x[:-1]))
    pair = array("I", [x[idx], x[idx + 1]])

    out = cc.substitute_pairs(x[:], pair, sub_value)

    assert len(out) < len(x)
    assert max(out) > max(x)


@given(generate_sequence())
def test_truncation(inputs):
    x, order = inputs

    etc_vf = cetc.compute(x, order, verbose=True, truncate=False)["ETC1D"]
    etc_vt = cetc.compute(x, order, verbose=True, truncate=True)["ETC1D"]
    etc_cf = cetc.compute(x, order, verbose=False, truncate=False)["ETC1D"]
    etc_ct = cetc.compute(x, order, verbose=False, truncate=True)["ETC1D"]

    assert etc_vf == etc_vt == etc_cf == etc_ct


def test_compute_save(tmp_path):

    x = array("I", [2, 4] * 100)
    file = tmp_path / "test.csv"

    etc_vf = cetc.compute_save(x, file, order=2, truncate=False)
    assert isinstance(etc_vf, dict)

    etc_vt = cetc.compute_save(x, file, order=2, truncate=True)
    assert isinstance(etc_vt, dict)
