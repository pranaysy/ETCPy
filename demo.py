#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: pranay
"""

# Import statements
from ETC.ETC import ETC, nETC
from ETC.utils import generate

# Generate synthetic data from the discrete uniform distribution
synthetic_data = generate(size=100, partitions=3)

# Compute Effort To Compress using Non-Sequential Recursive Pair Substitution
ETC(synthetic_data)

# Compute normalized version of ETC
nETC(synthetic_data)