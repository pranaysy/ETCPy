# -*- coding: utf-8 -*-

__import__("pkg_resources").declare_namespace(__name__)
import ETC.utils
from ETC.helper import parallel_1D, parallel_2D
from ETC.helper.IO import read, save
from ETC.NSRWS1D.compute_etc import compute as compute_1D
from ETC.NSRWS1D.compute_etc import compute_save as compute_save_1D
from ETC.NSRWS1D.compute_one_step import one_step as one_step_1D

from ETC.NSRWS2D.compute_etc import compute as compute_2D
from ETC.NSRWS2D.compute_etc import compute_save as compute_save_2D
from ETC.NSRWS2D.compute_one_step import one_step as one_step_2D

from ETC.complexity.compute_CCC import compute as compute_CCC

from ETC.utils import generate, partition
