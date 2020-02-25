#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: pranay
"""

from .NSRPS import one_step
from .utils import warmup, check_equality

warmup()


def ETC(x):

    etc = 0

    while not check_equality(x):

        x = one_step(x)
        etc += 1

    print(f">> Effort To Compress (ETC) = {etc}")
    return etc


def nETC(x, etc=None):

    if etc:
        etc *= 1 / len(x)
        print(f">> Normalized Effort To Compress (nETC) = {etc}")
        return etc
    else:
        etc = ETC(x) / len(x)
        print(f">> Normalized Effort To Compress (nETC) = {etc}")
        return etc
