#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: pranay
"""

# Import statements
import ETC

# Generate synthetic data from the discrete uniform distribution
synthetic_data = ETC.utils.generate(size=100, partitions=4)

# Compute Effort To Compress using Non-Sequential Recursive Pair Substitution
out = ETC.compute(synthetic_data, verbose=True)

# from multiprocessing import Pool
# from pathlib import Path

# folder = Path("/home/pranay/Projects/GenomeComplexity/data/GISAID/sequence")
# files = folder.glob("*.txt")

# def comp(filepath):
#     seq = ETC.IO.read(filepath)
#     fname = filepath.with_name(filepath.stem+'_etc_ord2.csv')
#     return filepath.stem, ETC.compute_save(seq, fname, order=2)

# def pcompute(z):
#     pool = Pool()
#     x = pool.map_async(comp, z)
#     pool.close()
#     pool.join()

#     return x

# if __name__ == '__main__':
#     a = pcompute(files)

