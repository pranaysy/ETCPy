#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""
import ETC
from ETC.Causal.pairs import ETC_causality, LZ_causality
from multiprocessing import Pool

def kernel_seq(inputs):

    idx, seqs = inputs
    seq_x, seq_y = seqs

    out = {"index":idx}

    out.update(ETC_causality(seq_x, seq_y))
    print(".", end="")
    return out

def kernel_seq_LZ(inputs):

    idx, seqs = inputs
    seq_x, seq_y = seqs

    out = {"index":idx}

    out.update(LZ_causality(seq_x, seq_y))
    print(".", end="")
    return out

def kernel_file(idx, seq_filenames):

    seq_x_filename, seq_y_filename = seq_filenames
    out = {"index":idx, "seq_x":seq_x_filename.name, "seq_y":seq_y_filename.name}

    seq_x = ETC.recode_lexical(ETC.read(seq_x_filename))
    seq_y = ETC.recode_lexical(ETC.read(seq_y_filename))

    out.update(ETC_causality(seq_x, seq_y))
    print(".", end="")
    return out

def parallelized(pairs, kernel):
    """
    This function operates concurrently on a collection of sequences. Loads
    each sequence and computes ETC.

    CAUTION: main module is unguarded, do not run these functions as is,
        particularly on Windows.

    Parameters
    ----------
    iterable : list/tuple/generator
        Collection of integer sequences.

    Returns
    -------
    list of dict elements
        Each dictionary element contains index, length of sequence & ETC.

    """
    # Initialize pool of parallel workers
    pool = Pool()

    # Map-execute function across sequences
    out = pool.map_async(kernel, enumerate(pairs))

    # Graceful exit
    pool.close()
    pool.join()

    # Return collected results
    return out.get()


def megapar(values, kernel):
    """
    This function operates concurrently on a collection of sequences. Loads
    each sequence and computes ETC.

    CAUTION: main module is unguarded, do not run these functions as is,
        particularly on Windows.

    Parameters
    ----------
    iterable : list/tuple/generator
        Collection of integer sequences.

    Returns
    -------
    list of dict elements
        Each dictionary element contains index, length of sequence & ETC.

    """
    # Initialize pool of parallel workers
    pool = Pool()

    # Map-execute function across sequences
    out = pool.map(kernel, enumerate(values))

    # Graceful exit
    pool.close()
    pool.join()

    # Return collected results
    return out