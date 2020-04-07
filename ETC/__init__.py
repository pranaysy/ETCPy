# -*- coding: utf-8 -*-

__import__("pkg_resources").declare_namespace(__name__)
import ETC.utils
from ETC.NSRWS import run_once_NSRWS
from ETC.ETC import compute, compute_save
from ETC.utils import generate, partition
from ETC.IO import read, save