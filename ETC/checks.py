#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 19:01:38 2020

@author: pranay
"""


def check_inputs(x):

    if isinstance(x, (list, tuple)):
        if len(x) > 1 and all(isinstance(a, (int, float)) for a in x):
            print(">> Input is a valid collection with numeric elements  ...")
            return True
        else:
            print(">> Input has invalid type or number of elements ...")
            return False

    else:
        print(">> Input is not a tuple or list object ...")
        return False