#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""
from functools import partial
from itertools import islice

# Import functions from standard library modules
from multiprocessing import Pool

# Import local modules
import ETC
from ETC.seq.process import entropy

# from ETC.seq.process import entropy
import numpy as np

# Function definitions
def _compute_single_file(filepath, order=2):
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
    seq = ETC.seq.IO.read(filepath)
    seq = ETC.seq.recode.recode_lexical(seq)

    # Filename for writing output of ETC computation
    fname = filepath.with_name(filepath.stem + f"_ETC_order{order}.csv")

    # Prepare output dictionary
    out = {"file": filepath.name, "length": len(seq), "entropy": entropy(seq)}

    # Compute ETC, write to file and update output dictionary
    out.update(ETC.NSRWS.x1D.etc.compute_save(seq, fname, order=order, truncate=True))

    return out


def pcompute_files(filelist, order=2):
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
    func = partial(_compute_single_file, order=order)
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
    out = {"index": seq[0], "length": len(seq[1]), "entropy": entropy(seq[1])}

    # Compute ETC and update output dictionary
    out.update(ETC.compute_1D(seq[1], order=2, verbose=False, truncate=True))

    return out


# TEMPORARY MOD FOR ASSEMBLY-FREE
# def _compute_single_seq(seq):
#     """
#     This function operates on a single sequence and computes ETC.

#     Parameters
#     ----------
#     seq : tuple of 2 elements
#         1st element is index for tracking.
#         2nd element is a sequence of integers used for ETC computation.
#         Output of enumerate.

#     Returns
#     -------
#     out : dict
#         index of sequence, length of sequence and ETC estimate.

#     """
#     data = ETC.seq.recode.recode_lexical(seq[1])

#     # Prepare output dictionary
#     out = {"item": seq[0], "length": len(data)}

#     # Compute ETC and update output dictionary
#     out.update(ETC.compute_1D(data, order=2, verbose=False, truncate=True))
#     out.update({"Entropy": entropy(data, legacy=False)})

#     return out


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


def pcompute_numpy(nparr):
    """
    This function operates concurrently row-wise on a 2D NumPy array. Loads
    each sequence and computes ETC.

    CAUTION: main module is unguarded, do not run these functions as is,
        particularly on Windows.

    Parameters
    ----------
    nparr : numpy array, int, 2D
        Sequence present as column, each row representing a different sequence

    Returns
    -------
    list of dict elements
        Each dictionary element contains index, length of sequence & ETC.

    """
    assert (
        isinstance(nparr, np.ndarray) and nparr.ndim == 2 and nparr.dtype == np.uint32
    ), ">ERROR: Input must be 2D NumPy array of 32-bit unsigned integers (np.uint32)"
    # Initialize pool of parallel workers
    pool = Pool()

    # Map-execute function across sequences
    out = pool.map_async(
        _compute_single_seq, enumerate([nparr[idx] for idx in range(nparr.shape[0])])
    )

    # Graceful exit
    pool.close()
    pool.join()

    # Return collected results
    return out.get()
