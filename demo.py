#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is a demo script for showcasing this package's functionality in brief.

@author: Pranay S. Yadav
"""

# Import call
import ETC


# IO & SEQUENCE MANAGEMENT
# ------------------------
# Read data to a list
text = ETC.read(filepath="somefile.txt", delimiter=",")

# Recode data to integers in lexicographic order
ETC.recode_lexical("bbacdbedf", case_sensitive=False)

# Check validity of input and automatically cast to the right form if valid
ETC.cast(text)

# Partition real-valued data to integer-valued discrete data
ETC.partition([0.1, 0.34, 0.68, -1.9, 25.3], n_bins=2)

# Generate synthetic data from the discrete uniform distribution
ETC.generate(size=1000, partitions=4)

# Compute Shannon Entropy for a sequence
ETC.entropy(seq=[1,2,1,1,1,2,1])


# 1D ETC ESTIMATION FOR A SINGLE SEQUENCE
# ---------------------------------------
# Generate a random discrete symbolic sequence
seq = ETC.generate(size=1000, partitions=2)

# Compute Effort To Compress using Non-Sequential Recursive Pair Substitution
out = ETC.compute_1D(seq, order=2, verbose=True, truncate=True)

# The result is a dictionary of substitutions plus other features
print(out)

# It can be saved to CSV using a convenience function:
ETC.save(out, filename="ETC_results.csv")


# PARALLELIZED 1D ETC ESTIMATION FOR CHUNKS OF A SINGLE SEQUENCE
# --------------------------------------------------------------
# Generate a long sequence
seq = ETC.generate(size=20000, partitions=2)

# Compute ETC on overlapping strides of 1000 elements offsetted by 100, in parallel
if __name__ == "__main__":
    outp = ETC.pcompute_single(seq, size=1000, offset=100)

# Compute ETC on non-overlapping strides of 1000 elements (set offset = size), in parallel
if __name__ == "__main__":
    outp = ETC.pcompute_single(seq, size=1000, offset=1000)


# PARALLELIZED 1D ETC ESTIMATION  FOR MULTIPLE SEQUENCES IN PARALLEL
# ------------------------------------------------------------------
# Generate 10 long sequences
seqs = (ETC.generate(20000, 4) for _ in range(10))

# Compute ETC estimates for each sequence
if __name__ == "__main__":
    outps = ETC.pcompute_multiple_seq(seqs)


# 2D ETC ESTIMATION FOR A PAIR OF SEQUENCES
# -----------------------------------------
# Generate two random sequences
seq_x = ETC.generate(size=1000, partitions=2)
seq_y = ETC.generate(size=1000, partitions=2)

# Compute Effort To Compress using Non-Sequential Recursive Pair Substitution
out = ETC.compute_2D(seq_x, seq_y, order=2, verbose=True, truncate=True)


# CCC ESTIMATION
# --------------
# Import call for CCC sub-package
from ETC import CCC

# Compute CCC for the above two sequences
ccc_est = CCC.compute(
    seq_x, seq_y, LEN_past=150, ADD_meas=15, STEP_size=20, n_partitions=False
)

# See docstrings for more information on CCC estimation
# ?CCC.compute