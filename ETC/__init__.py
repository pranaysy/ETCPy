# -*- coding: utf-8 -*-

__import__("pkg_resources").declare_namespace(__name__)

from ETC.seq.IO import read, save
from ETC.seq.process import generate, entropy
from ETC.seq.recode import cast, recode_lexical, partition, partition_numpy
from ETC.seq import check

from ETC.NSRWS.x1D.etc import compute as compute_1D

# from ETC.NSRWS.x1D.etc import compute_save as compute_save_1D
# from ETC.NSRWS.x1D.onestep import onestep as onestep_1D

from ETC.NSRWS.x2D.etc import compute as compute_2D

# from ETC.NSRWS.x2D.etc import compute_save as compute_save_2D
# from ETC.NSRWS.x2D.onestep import onestep as onestep_2D

# from ETC.CCC.compute_CCC import compute as compute_CCC


from ETC.NSRWS.x1D.parallel import (
    pcompute_multiple_seq,
    pcompute_single,
    pcompute_files,
    pcompute_numpy,
)

from ETC.LZ76.lzc import compute_complexity as LZC
from ETC.CCMC.pairs import CCM_causality
from ETC.CCMC.pairs_parallel import parallelized as CCM_causality_parallel
from ETC.CCMC.pairs_parallel import get_rowpairs
