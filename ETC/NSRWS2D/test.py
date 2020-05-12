#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""


from array import array
from random import choices, seed

from ETC.NSRWS2D.compute_core import check_equality as CE
from ETC.NSRWS2D.compute_etc import compute
from ETC.NSRWS2D.compute_one_step import one_step

seed(1)
size = 20
x = array("I", choices([1, 2], k=size))
y = array("I", choices([1, 2], k=size))

# x = array("I", [1]*10+[2])
# y = array("I", [1]*10+[2])

# out = one_step(x, y, 2, 1, 0)
# etc_x = cvtx(x)
# etc_y = cvtx(y)

trunc = compute(seq_x, seq_y, 2, 1, 0)
