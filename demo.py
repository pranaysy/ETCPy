#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is a demo script for showcasing this package's functionality in brief.

@author: Pranay S. Yadav
"""

# Import statements
import ETC

# 1D ETC ESTIMATION FOR A SINGLE SEQUENCE
# ---------------------------------------
# Generate synthetic data from the discrete uniform distribution
synthetic = ETC.generate(size=1000, partitions=4)

# Compute Effort To Compress using Non-Sequential Recursive Pair Substitution
out = ETC.compute_1D(synthetic, order=2, verbose=True, truncate=True)


# PARALLELIZED 1D ETC ESTIMATION FOR CHUNKS OF A SINGLE SEQUENCE
# --------------------------------------------------------------
# Generate a long sequence
seq = ETC.generate(10000, 4)

# Compute strides of 1000 elements offsetted by 100, in parallel
if __name__ == "__main__":
    outp = ETC.parallel_1D.pcompute_single(seq, size=1000, offset=100)

# Compute non-overlapping strides of 1000 elements (set offset = size)
if __name__ == "__main__":
    outp = ETC.parallel_1D.pcompute_single(seq, size=1000, offset=1000)

# PARALLELIZED 1D ETC ESTIMATION  FOR MULTIPLE SEQUENCES IN PARALLEL
# ------------------------------------------------------------------
# Generate 10 long sequences
seqs = (ETC.generate(10000, 4) for _ in range(10))

# Compute ETC estimates for each sequence
if __name__ == "__main__":
    outps = ETC.parallel_1D.pcompute_multiple_seq(seqs)


# 2D ETC ESTIMATION FOR A PAIR OF SEQUENCES
# -----------------------------------------
# Generate synthetic data from the discrete uniform distribution
seq_x = ETC.generate(size=1000, partitions=4)
seq_y = ETC.generate(size=1000, partitions=4)

# Compute Effort To Compress using Non-Sequential Recursive Pair Substitution
out = ETC.compute_2D(seq_x, seq_y, order=2, verbose=True, truncate=True)

# CCC ESTIMATION
# --------------
# Past window size
LEN_past = 150

# Current window size
ADD_meas = 15

# Jump step size
STEP_size = 20

ccc_est = ETC.compute_CCC(seq_x, seq_y, LEN_past, ADD_meas, STEP_size, partitions=False)

