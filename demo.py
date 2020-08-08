#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is a demo script for showcasing this package's functionality in brief.

@author: Pranay S. Yadav
"""

# Import call
import ETC

# ------------------------
# IO & SEQUENCE MANAGEMENT
# ------------------------
# Read data to a list
text = ETC.read(filepath="somefile.txt", delimiter=",") # Pick any file

# Check validity of input and automatically cast to the right form if valid
ETC.cast(text)

# Recode data to integers in lexicographic order
ETC.recode_lexical("bbacdbedf", case_sensitive=False)

# Partition real-valued data to integer-valued discrete data
ETC.partition([0.1, 0.34, 0.68, -1.9, 25.3], n_bins=2)

# Generate synthetic data from the discrete uniform distribution
ETC.generate(size=1000, partitions=4)

# Reproducibility of random generation can be controlled by passing the same seed value
ETC.generate(size=1000, partitions=4, seed=101)

# Compute Shannon Entropy for a sequence
ETC.entropy(seq=[1, 2, 1, 1, 1, 2, 1])


# ---------------------------------------
# 1D ETC ESTIMATION FOR A SINGLE SEQUENCE
# ---------------------------------------
# Generate a random discrete symbolic sequence
seq = ETC.generate(size=1000, partitions=2, seed=31)

# Simplest way to run
out = ETC.compute_1D(seq)

# The result is a dict of 2 key-value pairs: the raw and normalized ETC estimates
print(out)

# Get whichever is needed by using their respective keys
print(out.get('ETC1D'))
# [Out]: 225

print(out.get('NETC1D'))
# [Out]: 0.22522522522522523

# The normalization is done over one less than the overall length
print(out.get('ETC1D') / (len(seq) - 1))
# [Out]: 0.22522522522522523

# If more details about the trajectory are desired, set verbosity to True
out = ETC.compute_1D(seq, verbose=True)

# The result is now a dict of 3 elements: the 2 ETC estimates and the Trajectory
print(out.get('Trajectory')) # List of dicts - one dict for each step

# The default behavior is to truncate the iteration process until the sequence gets
# saturated to have all unique pairs occurring just once. This speeds up computation as
# the remaining steps don't need to be computed and ETC reduces to an analytic expression.
# However, the substitution table or features may be of interest and this truncation can
# then be turned off so that the iteration continues till entropy of 0 is reached:
out = ETC.compute_1D(seq, verbose=True, truncate=False)

print(out.get('Trajectory')) # Last step has length 1 and entropy 0

# This Trajectory can be saved to CSV for later use through a convenience function:
ETC.save(out.get('Trajectory'), filename="ETC_results.csv")

# -------------------------------------------------------------------------------------#
# Additionally, instead of pair-substitution (NSRPS), a window of any size may be
# substituted using the order switch, for example substitute triplets:
out = ETC.compute_1D(seq, order=3, verbose=True, truncate=False)

print(out.get("Trajectory"))

# The default function call ETC.compute_1D(seq) is the same as:
# ETC.compute_1D(seq, order=2, verbose=False, truncate=True)

# --------------------------------------------------------------
# PARALLELIZED 1D ETC ESTIMATION FOR CHUNKS OF A SINGLE SEQUENCE
# --------------------------------------------------------------
# Generate a long sequence
seq = ETC.generate(size=20000, partitions=2)

# Compute ETC on overlapping chunks of 1000 elements offsetted by 100, in parallel
if __name__ == "__main__":
    outp = ETC.pcompute_single(seq, size=1000, offset=100)

# The output is a list of dictionaries with estimates, one dict for each ordered chunk
print(outp)

# Compute ETC on non-overlapping chunks of 1000 elements (set offset = size), in parallel
if __name__ == "__main__":
    outp = ETC.pcompute_single(seq, size=1000, offset=1000)

# Similarly,
print(outp)

# ------------------------------------------------------------------
# PARALLELIZED 1D ETC ESTIMATION  FOR MULTIPLE SEQUENCES IN PARALLEL
# ------------------------------------------------------------------
# Generate 10 long sequences
seqs = [ETC.generate(10000, 2) for _ in range(10)]

# Compute ETC estimates for each sequence in parallel
if __name__ == "__main__":
    outp = ETC.pcompute_multiple_seq(seqs)

print(outp)


# --------------------------------
# WORKS WITH NUMPY OUT OF THE BOX!
# --------------------------------
# Generate a random discrete symbolic sequence and compute 1D ETC on it
import numpy as np
np.random.seed(10)
seq = np.random.randint(1, 3, size=5000)
out = ETC.compute_1D(seq)

print(out)
# {'ETC1D': 884, 'NETC1D': 0.17683536707341468}

# Parallelized ETC estimation - row-wise for 2D numpy arrays
seq = np.random.normal(1, 3, size=[10,5000]) # Each row is a distinct sequence
seq = ETC.partition_numpy(nparr=seq, n_bins=2)
out = ETC.pcompute_numpy(nparr=seq)

print(out)
# One estimate per row

# -----------------------------------------
# 2D ETC ESTIMATION FOR A PAIR OF SEQUENCES
# -----------------------------------------
# Generate two random sequences
seq_x = ETC.generate(size=1000, partitions=2, seed=17)
seq_y = ETC.generate(size=1000, partitions=2, seed=19)

# Compute Effort To Compress using Non-Sequential Recursive Pair Substitution
out = ETC.compute_2D(seq_x, seq_y, order=2, verbose=True, truncate=False)

# View estimates
print(out.get('ETC2D'))
print(out.get('NETC2D'))

# View trajectory
print(out.get('Trajectory'))

# -----------------------------------------
# CAUSALITY TESTING USING THE CCC FRAMEWORK
# -----------------------------------------
# Import call for CCC sub-package
from ETC import CCC

# Compute CCC for the above two sequences
ccc_est = CCC.compute(
    seq_x, seq_y, LEN_past=150, ADD_meas=15, STEP_size=20, n_partitions=False
)
# [Out]: CCC for seq_y -> seq_x = -0.00301035257856264

# See docstrings for more information on CCC estimation
# ?CCC.compute

# Simulate a pair of coupled first-order AR processes
ar = CCC.coupled_AR(length=10000, a=0.8, b=0.9, c=0.8, e=0.01, burn=1000, seed=1)
# ar is a dictionary of two key-value pairs with the following keys:
#   "dependent" and "independent", each with their respective values in float arrays
# ?CCC.coupled_AR for more information on sampling from AR processes

# Estimate CCC for the direction independent -> dependent with binning
ccc_ar = CCC.compute(
    seq_x=ar["dependent"],
    seq_y=ar["independent"],
    LEN_past=150,
    ADD_meas=15,
    STEP_size=20,
    n_partitions=2,
)
# [Out]: CCC for seq_y -> seq_x = 0.005755172746030292

# And for the opposite direction
ccc_ar = CCC.compute(
    seq_x=ar["independent"],
    seq_y=ar["dependent"],
    LEN_past=150,
    ADD_meas=15,
    STEP_size=20,
    n_partitions=2,
)
# [Out]: CCC for seq_y -> seq_x = 0.0002971309733327245
