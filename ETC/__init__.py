# -*- coding: utf-8 -*-

__import__("pkg_resources").declare_namespace(__name__)

import logging

logger = logging.getLogger(__name__)

if logger.hasHandlers():
    logger.handlers.clear()

c_handler = logging.StreamHandler()
c_handler.setLevel("DEBUG")
c_format = logging.Formatter(
    fmt="{asctime} {levelname} {module}:{lineno:03d} {message} @{funcName}",
    datefmt="%d-%b-%y %H:%M:%S",
    style="{",
)
c_handler.setFormatter(c_format)

logger.addHandler(c_handler)
logger.setLevel("DEBUG")


logger.info("And it begins")

from ETC.utils.IO import read, save
from ETC.utils.process import generate, entropy
from ETC.utils.recode import cast, recode_lexical, partition, partition_numpy
from ETC.utils import check

from ETC.NSRWS.x1D.etc import compute as compute_1D

# from ETC.NSRWS.x1D.etc import compute_save as compute_save_1D
# from ETC.NSRWS.x1D.onestep import onestep as onestep_1D

from ETC.NSRWS.x2D.etc import compute as compute_2D

# from ETC.NSRWS.x2D.etc import compute_save as compute_save_2D
# from ETC.NSRWS.x2D.onestep import onestep as onestep_2D

# from ETC.CCC.compute_CCC import compute as compute_CCC
from ETC.LZ76.lzc import compute_complexity as LZC

logger.info("Core module import successful")

from ETC.NSRWS.x1D.parallel import (
    pcompute_multiple_seq,
    pcompute_single,
    pcompute_files,
    pcompute_numpy,
)

from ETC.CCMC.pairs import CCM_causality
from ETC.CCMC.pairs_parallel import parallelized as CCM_causality_parallel
from ETC.CCMC.pairs_parallel import get_rowpairs
