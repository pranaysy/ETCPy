#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""


from pathlib import Path

import ETC

cov2_file = Path(
    "/home/pranay/Projects/GenomeComplexity/data/SARS-CoV-1/refseq/NC_004718.3_RefSeq.txt"
)
cov1_file = Path(
    "/home/pranay/Projects/GenomeComplexity/data/SARS-CoV-2_old/refseq/NC_045512.2_RefSeq.txt"
)

cov2_bin4 = ETC.read(cov2_file)
cov1_bin4 = ETC.read(cov1_file)

out_bin4 = ETC.parallel.pcompute_multiple_seq((cov1_bin4, cov2_bin4))

cov1cov2_bin4 = cov1_bin4 + cov2_bin4
cov2cov1_bin4 = cov2_bin4 + cov1_bin4

out_bin4 = ETC.parallel.pcompute_multiple_seq(
    (cov1_bin4, cov2_bin4, cov1cov2_bin4, cov2cov1_bin4)
)

cov1_netc_bin4 = out_bin4[0]["ETC"] / (out_bin4[0]["length"] - 1)
cov2_netc_bin4 = out_bin4[1]["ETC"] / (out_bin4[1]["length"] - 1)
cov1cov2_netc_bin4 = out_bin4[2]["ETC"] / (out_bin4[2]["length"] - 1)
cov2cov1_netc_bin4 = out_bin4[3]["ETC"] / (out_bin4[3]["length"] - 1)

d_netc_bin4 = cov1cov2_netc_bin4 - cov1_netc_bin4 + cov2cov1_netc_bin4 - cov2_netc_bin4

cov2_bin2 = ETC.helper.IO.read_dna(cov2_file)
cov1_bin2 = ETC.helper.IO.read_dna(cov1_file)

cov1cov2_bin2 = cov1_bin2 + cov2_bin2
cov2cov1_bin2 = cov2_bin2 + cov1_bin2

out_bin2 = ETC.parallel.pcompute_multiple_seq(
    (cov1_bin2, cov2_bin2, cov1cov2_bin2, cov2cov1_bin2)
)

cov1_netc_bin2 = out_bin2[0]["ETC"] / (out_bin2[0]["length"] - 1)
cov2_netc_bin2 = out_bin2[1]["ETC"] / (out_bin2[1]["length"] - 1)
cov1cov2_netc_bin2 = out_bin2[2]["ETC"] / (out_bin2[2]["length"] - 1)
cov2cov1_netc_bin2 = out_bin2[3]["ETC"] / (out_bin2[3]["length"] - 1)

d_netc_bin2 = 0.5 * (
    cov1cov2_netc_bin2 - cov1_netc_bin2 + cov2cov1_netc_bin2 - cov2_netc_bin2
)

print(f"> Distance metric for ETC on 2 bins: {d_netc_bin2}")
print(f"> Distance metric for ETC on 4 bins: {d_netc_bin4}")

# %%
k = 100
cov1_bin2 = cov1_bin2[:k]
cov2_bin2 = cov2_bin2[:k]
cov1cov2_bin2 = cov1_bin2 + cov2_bin2
cov2cov1_bin2 = cov2_bin2 + cov1_bin2

out_bin2 = ETC.parallel.pcompute_multiple_seq(
    (cov1_bin2, cov2_bin2, cov1cov2_bin2, cov2cov1_bin2)
)

cov1_netc_bin2 = out_bin2[0]["ETC"] / (out_bin2[0]["length"] - 1)
cov2_netc_bin2 = out_bin2[1]["ETC"] / (out_bin2[1]["length"] - 1)
cov1cov2_netc_bin2 = out_bin2[2]["ETC"] / (out_bin2[2]["length"] - 1)
cov2cov1_netc_bin2 = out_bin2[3]["ETC"] / (out_bin2[3]["length"] - 1)

d_netc_bin2 = 0.5 * (
    cov1cov2_netc_bin2 - cov1_netc_bin2 + cov2cov1_netc_bin2 - cov2_netc_bin2
)

print(f"> Distance metric for ETC on 2 bins: {d_netc_bin2}")
print(f"> Distance metric for ETC on 4 bins: {d_netc_bin4}")

# %%
def runner(seq):
    out = ETC.compute(seq, 2, 0, 1)
    return out["ETC"] / (len(seq) - 1)


def test():
    a = ETC.generate(100, 2)
    b = ETC.generate(100, 2)
    ab = a + b
    ba = b + a
    d = 0.5 * (runner(ab) - runner(a) + runner(ba) - runner(b))
    print(d)
    if d > 0:
        return True
    else:
        return False
