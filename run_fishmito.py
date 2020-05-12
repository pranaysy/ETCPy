#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""
import itertools as it
from pathlib import Path

import ETC

filelist = ETC.helper.IO.populate_files(
    "/home/pranay/Projects/GenomeComplexity/data/AFProject/sequence/assembled-fish_mito"
)

for file1, file2 in it.combinations(filelist, 2):
    text1 = Path(file1).read_text()
    text2 = Path(file2).read_text()
    text12 = text1 + text2
    file12name = "cat_" + file1.stem + "_" + file2.stem + ".txt"
    file1.with_name(file12name).write_text(text12)
    text21 = text2 + text1
    file21name = "cat_" + file2.stem + "_" + file1.stem + ".txt"
    file2.with_name(file21name).write_text(text21)
    print(f"> Done: {file1.stem} & {file2.stem}")

filelist = ETC.helper.IO.populate_files(
    "/home/pranay/Projects/GenomeComplexity/data/AFProject/sequence/assembled-fish_mito"
)
filelist = tuple(filelist)

# if __name__ == '__main__':
#     results = ETC.parallel.pcompute_files(filelist, order=2)

for n_order in [2, 3, 4]:
    results = ETC.parallel.pcompute_files(filelist, order=n_order)
    ETC.helper.IO.save(results, f"ETC_estimates{n_order}.csv")

# %%

x = ETC.helper.IO.read(
    "/home/pranay/Projects/GenomeComplexity/data/AFProject/sequence/assembled-ecoli/D1Sd197_seq.txt"
)
