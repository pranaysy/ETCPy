# cython: language_level=3, boundscheck=False, wraparound=False, nonecheck=False, emit_code_comments=True, cdivision=True, embedsignature=True
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""
# Import stuff
from cpython cimport array, bool
cimport cython
import array

# Function for getting mask for pairs
cpdef array.array get_mask_pairs(const unsigned int[::1] x):
    """
    INPUT
    -----
    x : array.array
        Array object containing 32-bit integers.

    OUTPUT
    ------
    mask : array.array
        Array object containing 32-bit integers - 0s or 1s corresponding to values in
        x for which successive overlapping pairs occur.
    """
    # Get size of input
    cdef Py_ssize_t x_size = len(x)

    # Initialize a mask of Falses
    cdef array.array int_template = array.array('I', [])
    cdef array.array mask = array.clone(int_template, x_size-1, zero=True)
    cdef unsigned int[:] mask_view = mask

    # Initialize bounds for iteration
    cdef Py_ssize_t n = 0

    # Turn all values in mask to Trues
    for n in range(x_size-1):
        mask_view[n] += 1

    # Iterate over all values of input
    n = 0
    while n < x_size - 2:

        # If successive pairs match
        if x[n] == x[n+1] and x[n+1] == x[n+2]:

            # Mask out the second one
            mask_view[n+1] = 0

            # And slide over it
            n += 1

        # Increment while loop index
        n += 1

    return mask

# Function for substituting pairs
cpdef list substitute_pairs(unsigned int[::1] x, unsigned int[::1] pair, unsigned int value):
    """
    INPUT
    -----
    x : array.array
        Array object containing 32-bit unsigned integers.

    pair : array.array, length = 2
        Array object containing 2 32-bit unsigned integers.

    value : unsigned 32-bit int
        Value to substitute the first element of pair with

    OUTPUT
    ------
    out : list
        Array object containing 32-bit integers, with supplied pair replaced everywhere
        by the supplied value.
    """
    # Initialize looping variables and output list
    cdef Py_ssize_t n = 0
    cdef Py_ssize_t x_size = len(x)
    cdef list out = []

    # Loop over input and replace pair
    while n < x_size-1:

        # Check for match with supplied pair
        if x[n] == pair[0] and x[n+1] == pair[1]:

            # Replace first value with supplied value
            x[n] = value

            # Replace second value with 0
            x[n+1] = 0

            n += 1

        n += 1

    # Reset indexing variable
    n = 0

    # Loop over mutated input and append non-zero values to list
    for n in range(x_size):

        if x[n]:

            out.append(x[n])

    return out

# Function for checking whether all elements in input are identical
cpdef bint check_equality(const unsigned int[::1] x):
    """
    INPUT
    -----
    x : array.array
        Array object containing 32-bit unsigned integers.


    OUTPUT
    ------
    bool
        True if all elements are identical
    """
    # Intialize loop bounds
    cdef Py_ssize_t n
    cdef Py_ssize_t x_size = len(x)

    # Iterate over values from input
    for n in range(x_size):

        # Short-circuit the loop: check for any element that doesn't equal the first
        if x[0] != x[n]:
            return False

    return True

# Function for getting mask for windows of any length
cpdef array.array get_mask_windows(const unsigned int[::1] x, unsigned int order):
    """
    INPUT
    -----
    x : array.array
        Array object containing 32-bit integers.

    order: unsigned 32-bit int
        Length of the window to slide across input

    OUTPUT
    ------
    mask : array.array
        Array object containing 32-bit integers - 0s or 1s corresponding to values in
        x for which successive overlapping pairs occur.
    """
    # Get size of input
    cdef Py_ssize_t x_size = len(x)

    # Initialize a mask of Falses
    cdef array.array int_template = array.array('I', [])
    cdef array.array mask = array.clone(int_template, x_size - (order-1), zero=True)
    cdef unsigned int[:] mask_view = mask

     # Initialize variable for iteration
    cdef Py_ssize_t n = 0

    # Turn all values in mask to Trues
    for n in range(x_size - order + 1):
        mask_view[n] += 1

    # Initialize variables for iteration across input
    cdef Py_ssize_t k = 0 # Outer loop
    cdef Py_ssize_t m = 0 # Inner loop

    # Tracking variable for counting matching elements in pairwise window comparison
    cdef unsigned int track = 0

    # Iterate over input values except the last 'order' values [Outermost master loop]
    n = 0
    for n in range(x_size - (order-1)):

        # proceed only if mask is True for current element
        if mask_view[n]:

            # Outer loop for sliding the 'next' window by unit step (current vs next)
            for k in range(1,order):  # Start from 1 - begin comparing from next window

                # Inner loop for comparing elements in current and next windows
                for m in range(order):

                    if n+m+k >= x_size:
                        return mask

                    # If elements match, increment tracker
                    if x[n+m] == x[n+m+k]:
                        track += 1

                    # Else stop iteration over this comparison of windows
                    # else:
                    #     break

                # Trick: preserve mask only if track doesn't equal order
                # If track == order, short-circuit eval takes precedence, returning 0
                mask_view[n+k] = track!=order and mask_view[n+k]

                # Reset tracker
                track = 0

    return mask

# Function for substituting windows of any length
cpdef list substitute_windows(unsigned int[::1] x, unsigned int order, unsigned int[::1] window, unsigned int value):
    """
    INPUT
    -----
    x : array.array
        Array object containing 32-bit unsigned integers.

    order: unsigned 32-bit int
        Length of the window to slide across input

    window : array.array, length = 2
        Array object containing 2 32-bit unsigned integers.

    value : unsigned 32-bit int
        Value to substitute the first element of pair with

    OUTPUT
    ------
    out : list
        Array object containing 32-bit integers, with supplied pair replaced everywhere
        by the supplied value.
    """
    # Initialize looping variables and output list
    cdef Py_ssize_t n = 0 # Outer loop
    cdef Py_ssize_t m = 0 # Inner loop
    cdef Py_ssize_t x_size = len(x)
    cdef list out = []

    # Tracking variable for counting matching elements in pairwise window comparison
    cdef unsigned int track = 0

    # Iterate over input values except one less than the last 'order' values
    # Logic: last window, say triplet must begin from 3rd-last index, leaving 2 values
    for n in range(x_size-order+1):

        # Slide window of given order and do element-wise comparison
        for m in range(order):

            # Reset tracker


            # Track comparison of input elements with window elements
            if x[n+m] == window[m]:
                track += 1
            # # If mismatch, break
            else:
                break

        # If all compared elements match for current window
        if track == order:

            # Replace the first element with provided value
            x[n] = value

            # Replace the remaining subsequent values with zeros
            for m in range(1, order):
                x[n+m] = 0

        # Reset tracker
        track = 0
    # Reset indexing variable
    n = 0

    # Loop over mutated input and append non-zero values to list
    for n in range(x_size):

        if x[n]:

            out.append(x[n])

    return out