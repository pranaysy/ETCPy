#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""

# %%
from ETC.helper.compute_markov_transition_probs import compute, _read_sequence, _compute_transition_probs, _generate_overlaps, sample_sequence
from random import choices
from pathlib import Path
import ETC

# %%
file = Path("/home/pranay/Projects/GenomeComplexity/data/AFProject/sequence/assembled-fish_mito/NC_009057_seq.txt")
seq = _read_sequence(file)
est = compute(file,1, compact=True, flatten=False)

# %%
tail = seq[-2:]
for n in range(10000):
    last = tail[-2:]
    probs = est.loc[last,:]
    new = ''.join(choices(probs.index, weights = probs.values, k=1))
    tail+=new

# %%
sample_sequence(seq, order=1, size=50, sampler_seed=64)


# %%
file1 = Path("/home/pranay/Projects/GenomeComplexity/data/AFProject/sequence/assembled-fish_mito/NC_009057_seq.txt")
file2 = Path("/home/pranay/Projects/GenomeComplexity/data/AFProject/sequence/assembled-fish_mito/NC_009058_seq.txt")