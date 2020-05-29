#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""
from array import array
from hypothesis import given
from hypothesis import strategies as st
from ETC.seq import recode
from collections import Counter
from string import ascii_lowercase

invalid_types = (
    st.fractions(),
    st.characters(),
    st.floats(),
    st.text(),
    st.complex_numbers(),
    st.integers(max_value=0),
    st.integers(min_value=2 ** 32),
)
valid_types = st.integers(min_value=1, max_value=2 ** 32 - 1)

def counts(x):
    return tuple(Counter(x).values())

@given(
    x=st.one_of(st.tuples(st.one_of(invalid_types)), st.lists(st.one_of(invalid_types)))
)
def test_cast_invalid(x):

    x = recode.cast(x)

    assert x is None


def test_cast_zeroes():

    x = recode.cast([0, 0, 0, 0])

    assert x is None


@given(x=st.one_of(st.tuples(valid_types), st.lists(valid_types, min_size=1)))
def test_cast_valid(x):

    x = recode.cast(x)

    assert isinstance(x, array) and x.typecode == "I"

@given(x=st.text(min_size=1, alphabet=list(ascii_lowercase)))
def test_all_recodes(x):

    x1 = recode.recode_lexical(x)
    x2 = recode.recode_alphabetical(x)
    x3 = recode.recode_randint(x)
    x4 = recode.recode_random(x)

    assert counts(x1) == counts(x2) == counts(x3) == counts(x4)
    assert len(set(x1)) == len(set(x2)) == len(set(x3)) == len(set(x4))