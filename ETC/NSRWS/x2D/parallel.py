#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""
from functools import partial
from itertools import islice
from collections import Counter
from random import choices

# Import functions from standard library modules
from multiprocessing import Pool

# Import local modules
import ETC
from ETC.seq.process import entropy
from ETC.helper.compute_markov_transition_probs import sample_sequence

# Function definitions
def _compute_two_files_truncated(files, order=2):
    """
    This function operates on a single file - reads sequence, computes ETC
    and writes to disk.

    Parameters
    ----------
    filepath : str or Path object
        Valid path to a file containing sequence.

    Returns
    -------
    out : dict
        filename, length of sequence and ETC estimate.

    """
    # Read file as a sequence
    filepath1, filepath2 = files
    seq1 = ETC.helper.IO.read(filepath1)
    seq2 = ETC.helper.IO.read(filepath2)

    if len(seq1) > len(seq2):
        seq1 = seq1[: len(seq2)]
    else:
        seq2 = seq2[: len(seq1)]

    # Prepare output dictionary
    out = {"seq1": filepath1.stem, "seq2": filepath2.stem, "length": len(seq1)}

    # Compute ETC, write to file and update output dictionary
    out.update(ETC.compute_2D(seq1, seq2, order=order, truncate=True, verbose=False))
    seq1etc = ETC.compute_1D(seq1, order=order, truncate=True, verbose=False)["ETC1D"]
    out.update({"ETC1D_seq1": seq1etc})

    seq2etc = ETC.compute_1D(seq2, order=order, truncate=True, verbose=False)["ETC1D"]
    out.update({"ETC1D_seq2": seq2etc})

    return out


# Function definitions
def _compute_two_files_markov(files, markov_order, order=2):
    """
    This function operates on a single file - reads sequence, computes ETC
    and writes to disk.

    Parameters
    ----------
    filepath : str or Path object
        Valid path to a file containing sequence.

    Returns
    -------
    out : dict
        filename, length of sequence and ETC estimate.

    """
    # Read file as a sequence
    filepath1, filepath2 = files
    seq1 = ETC.helper.IO.read(filepath1, recode=False)
    seq2 = ETC.helper.IO.read(filepath2, recode=False)

    lseq1 = len(seq1)
    lseq2 = len(seq2)

    if lseq1 > lseq2:
        diff = lseq1 - lseq2
        extra_tail = sample_sequence(
            seq2, order=markov_order, size=diff, sampler_seed=64
        )
        seq2 += extra_tail

    elif lseq1 < lseq2:
        diff = lseq2 - lseq1
        extra_tail = sample_sequence(
            seq1, order=markov_order, size=diff, sampler_seed=64
        )
        seq1 += extra_tail

    assert len(seq1) == len(seq2)

    seq1 = ETC.helper.IO.recode_to_int(seq1)
    seq2 = ETC.helper.IO.recode_to_int(seq2)

    # Filename for writing output of ETC computation
    # fname = filepath1.with_name(filepath1.stem + '_&_'+ filepath2.stem + f"_etc_order{order}_markov_order{markov_order}.csv")

    # Prepare output dictionary
    out = {"seq1": filepath1.stem, "seq2": filepath2.stem, "length": len(seq1)}

    # Compute ETC, write to file and update output dictionary
    out.update(ETC.compute_2D(seq1, seq2, order=order, truncate=True, verbose=False))
    seq1etc = ETC.compute_1D(seq1, order=order, truncate=True, verbose=False)["ETC1D"]
    out.update({"ETC1D_seq1": seq1etc})

    seq2etc = ETC.compute_1D(seq2, order=order, truncate=True, verbose=False)["ETC1D"]
    out.update({"ETC1D_seq2": seq2etc})

    return out


def pcompute_files_markov(filelist, markov_order, order=2):
    """
    This function operates concurrently on a list of files. Reads each as a
    sequence, computes ETC and writes output to disk.

    CAUTION: main module is unguarded, do not run these functions as is,
        particularly on Windows.

    Parameters
    ----------
    filelist : list/tuple/generator
        Collection of filenames of files containing sequence data.

    Returns
    -------
    list of dict elements
        Each dictionary element contains filename, length of sequence & ETC.

    """
    # Initialize pool of parallel workers
    pool = Pool()
    func = partial(_compute_two_files_markov, markov_order=markov_order, order=order)
    # Map-execute function across files
    out = pool.map_async(func, filelist)

    # Graceful exit
    pool.close()
    pool.join()

    # Return collected results
    return out.get()


def pcompute_files_truncated(filelist, order=2):
    """
    This function operates concurrently on a list of files. Reads each as a
    sequence, computes ETC and writes output to disk.

    CAUTION: main module is unguarded, do not run these functions as is,
        particularly on Windows.

    Parameters
    ----------
    filelist : list/tuple/generator
        Collection of filenames of files containing sequence data.

    Returns
    -------
    list of dict elements
        Each dictionary element contains filename, length of sequence & ETC.

    """
    # Initialize pool of parallel workers
    pool = Pool()
    func = partial(_compute_two_files_truncated, order=order)
    # Map-execute function across files
    out = pool.map_async(func, filelist)

    # Graceful exit
    pool.close()
    pool.join()

    # Return collected results
    return out.get()


def _compute_single_seq(seq):
    """
    This function operates on a single sequence and computes ETC.

    Parameters
    ----------
    seq : tuple of 2 elements
        1st element is index for tracking.
        2nd element is a sequence of integers used for ETC computation.
        Output of enumerate.

    Returns
    -------
    out : dict
        index of sequence, length of sequence and ETC estimate.

    """
    # Prepare output dictionary
    out = {"item": seq[0], "length": len(seq[1]), "entropy": entropy(seq[1])}

    # Compute ETC and update output dictionary
    out.update(ETC.compute(seq[1], order=2, verbose=False, truncate=True))

    return out


def pcompute_multiple_seq(iterable):
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
    out = pool.map_async(_compute_single_seq, enumerate(iterable))

    # Graceful exit
    pool.close()
    pool.join()

    # Return collected results
    return out.get()


def _overlapping_chunks(seq, size, offset=1):
    """
    This function takes an input sequence and produces chunks of chosen size.
    Offset can be used to control degree of overlap (or distance between chunks
    that don't overlap)

    Parameters
    ----------
    seq : tuple or list
        Sequence of integers.
    size : int
        Length of each produced chunk.
    offset : int, optional
        Number of elements to shift each chunk by. The default is 1.
        Setting this to any value less than size allows control of overlap.
        Setting this >= size produces non-overlapping chunks.

    Returns
    -------
    zip
        zip object that produces chunks of specified size, one at a time.

    """

    return zip(*(islice(seq, i, None, offset) for i in range(size)))


def _non_overlapping_chunks(seq, size):
    """
    This function takes an input sequence and produces chunks of chosen size
    that strictly do not overlap. This is a much faster implemetnation than
    _overlapping_chunks and should be preferred if running on very large seq.

    Parameters
    ----------
    seq : tuple or list
        Sequence of integers.
    size : int
        Length of each produced chunk.

    Returns
    -------
    zip
        zip object that produces chunks of specified size, one at a time.

    """

    return zip(*[iter(seq)] * size)


def pcompute_single(seq, size, offset=1):
    """
    This function operates concurrently on chunks of a given sequence. Gets
    each chunk and computes ETC one-by-one. Offset parameter controls degree of
    overlap (or non-overlap)

    CAUTION: main module is unguarded, do not run these functions as is,
        particularly on Windows.

    Parameters
    ----------
    seq : tuple or list
        Sequence of integers.
    size : int
        Length of each produced chunk.
    offset : int, optional
        Number of elements to shift each chunk by. The default is 1.
        Setting this to any value less than size allows control of overlap.
        Setting this >= size produces non-overlapping chunks.

    Returns
    -------
    list of dict elements
        Each dictionary element contains index, length of sequence & ETC.

    """
    # If offset equals size, get non-overlapping chunks of given size
    if offset == size:
        iterable = _non_overlapping_chunks(seq, size)

    # Else get overlapping chunks of given size and offset
    else:
        iterable = _overlapping_chunks(seq, size, offset)

    # Execute parallel computation over chunks
    return pcompute_multiple_seq(iterable)
