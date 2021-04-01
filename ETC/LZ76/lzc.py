#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""
from ETC.LZ76 import core
from ETC.utils.recode import cast
from ETC.utils.check import arraytype

# from ETC import logger
import logging

logger = logging.getLogger(__name__)


def compute_complexity(seq):

    logger.info("LZC called")
    # Coerce input to appropriate array type, if not possible throw a fit & exit
    if not arraytype(seq):
        seq = cast(seq)
        if seq is None:
            return None

    # Check whether all elements are equal, & exit if True (LZ76 of such inputs is 2)
    if core.check_equality(seq):
        print("> All elements in sequence are equal!")
        return 2

    # Else execute Cython function for computing LZ complexity
    return core.lzc_a(seq)
