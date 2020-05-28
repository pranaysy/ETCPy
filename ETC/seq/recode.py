#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""

from string import ascii_lowercase
from random import shuffle, choices
from array import array


def cast(seq):

    out = array("I", seq)
    return out


def recode_lexical(text):

    alphabets = sorted(set(text))
    replacer = dict((y, str(x + 1)) for x, y in enumerate(alphabets))
    text = cast(int(replacer[x]) for x in text)
    return text


def recode_alphabetical(text):

    replacer = dict((y, str(x + 1)) for x, y in enumerate(ascii_lowercase))
    text = cast(int(replacer[x]) for x in text.lower())
    return text


def recode_dna(text):

    replacer = {"A": "1", "G": "1", "C": "2", "T": "2"}
    text = cast(int(replacer[x]) for x in text.upper())
    return text


def recode_random(text):

    alphabets = list(set(text))
    shuffle(alphabets)
    replacer = dict((y, str(x + 1)) for x, y in enumerate(alphabets))
    text = cast(int(replacer[x]) for x in text)
    return text


def recode_randint(text):

    alphabets = list(set(text))
    numbers = choices(range(1, 2 ** 20), k=len(alphabets))
    replacer = dict(zip(alphabets, numbers))
    text = cast(int(replacer[x]) for x in text)
    return text
