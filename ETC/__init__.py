# -*- coding: utf-8 -*-

__import__("pkg_resources").declare_namespace(__name__)
import ETC.utils
from ETC.helper import parallel
from ETC.helper.IO import read, save
from ETC.NSRWS1D.compute_etc import compute, compute_save
from ETC.NSRWS1D.compute_one_step import one_step
from ETC.utils import generate, partition
