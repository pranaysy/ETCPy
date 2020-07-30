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

get1D = partial(ETC.compute_1D, order=2, verbose=False, truncate=True)


def _compute_distance(inputs):
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
    idx, seqs = inputs

    S1 = ETC.seq.recode.recode_lexical(seqs[0])
    S2 = ETC.seq.recode.recode_lexical(seqs[1])

    # Prepare output dictionary
    out = {"item": idx, "length_seq1": len(S1), "length_seq2": len(S2)}


    # Compute ETC and update output dictionary
    ETC1D_seq1 = get1D(S1)["ETC1D"]
    out.update({"ETC1D_seq1": ETC1D_seq1})

    ETC1D_seq2 = get1D(S2)["ETC1D"]
    out.update({"ETC1D_seq2": ETC1D_seq2})

    ETC1D_seq1seq2 = get1D(S1 + S2)["ETC1D"]
    out.update({"ETC1D_seq1seq2": ETC1D_seq1seq2})

    ETC1D_seq2seq1 = get1D(S2 + S1)["ETC1D"]
    out.update({"ETC1D_seq2seq1": ETC1D_seq2seq1})

    dETC = 0.5 * (ETC1D_seq1seq2 + ETC1D_seq2seq1 - ETC1D_seq1 - ETC1D_seq2)

    out.update({"distance": dETC})

    return out


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
    out = pool.map(_compute_distance, enumerate(iterable))

    # Graceful exit
    pool.close()
    pool.join()

    # Return collected results
    return out


def truncate(seq1, seq2):

    # Truncate the longer sequence
    if len(seq1) == len(seq2):
        return seq1, seq2

    if len(seq1) > len(seq2):
        seq1 = seq1[: len(seq2)]
    else:
        seq2 = seq2[: len(seq1)]

    return seq1, seq2


def pcompute_single(seq1, seq2, size, offset=1):
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

    seq1, seq2 = truncate(seq1, seq2)

    # If offset equals size, get non-overlapping chunks of given size
    if offset == size:
        iterable1 = _non_overlapping_chunks(seq1, size)
        iterable2 = _non_overlapping_chunks(seq2, size)

    # Else get overlapping chunks of given size and offset
    else:
        iterable1 = _overlapping_chunks(seq1, size, offset)
        iterable2 = _overlapping_chunks(seq2, size, offset)

    # Execute parallel computation over chunks
    return pcompute_multiple_seq(zip(iterable1, iterable2))
