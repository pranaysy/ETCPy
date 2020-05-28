# -*- coding: utf-8 -*-

__import__("pkg_resources").declare_namespace(__name__)


from ETC.seq.IO import read, save
from ETC.NSRWS.x1D.etc import compute as compute_1D

# from ETC.NSRWS.x1D.etc import compute_save as compute_save_1D
# from ETC.NSRWS.x1D.onestep import onestep as onestep_1D

from ETC.NSRWS.x2D.etc import compute as compute_2D

# from ETC.NSRWS.x2D.etc import compute_save as compute_save_2D
# from ETC.NSRWS.x2D.onestep import onestep as onestep_2D

# from ETC.CCC.compute_CCC import compute as compute_CCC

from ETC.seq.process import generate, partition
from ETC.seq.recode import cast, recode_lexical


from ETC.NSRWS.x1D.parallel import pcompute_multiple_seq
