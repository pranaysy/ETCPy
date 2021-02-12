# cython: language_level=3, boundscheck=False, wraparound=False, nonecheck=False, emit_code_comments=True, cdivision=True, embedsignature=True
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""
# Import stuff
# cimport cython

cpdef unsigned int lzc_a(const unsigned int[::1] intarray):
    """
    Lempel-Ziv (LZ76) complexity on 32-bit integer arrays
    """

    # Variables Initialization
    cdef Py_ssize_t arraylength = len(intarray)
    cdef unsigned int complexity = 1
    cdef Py_ssize_t prefix_len = 1
    cdef Py_ssize_t len_substring = 1
    cdef Py_ssize_t max_len_substring = 1
    cdef unsigned int pointer = 0

    # While we haven't decoded the full string we continue
    while prefix_len + len_substring <= arraylength:

        # Given a prefix length, find the largest substring
        if (
            intarray[pointer + len_substring - 1]
            == intarray[prefix_len + len_substring - 1]
        ):
            len_substring += 1  # increase the length of the substring
        else:

            max_len_substring = max(len_substring, max_len_substring)
            pointer += 1

            # all the pointers have been investigated, we pick the largest for the jump
            if pointer == prefix_len:

                # Increment complexity
                complexity += 1

                # Increase the prefix length by the maximum substring size found so far
                prefix_len += max_len_substring

                # Reset the variables
                pointer = 0
                max_len_substring = 1

            # reset the length of the substring
            len_substring = 1

    # Check final repetition if we were in the middle of a substring
    if len_substring != 1:
        complexity += 1

    return complexity

cpdef unsigned int lzc_b(const unsigned char[:] bytestring):
    """
    Lempel-Ziv (LZ76) complexity on bytestrings
    """

    # Variables Initialization
    cdef Py_ssize_t stringlength = len(bytestring)
    cdef unsigned int complexity = 1
    cdef Py_ssize_t prefix_len = 1
    cdef Py_ssize_t len_substring = 1
    cdef Py_ssize_t max_len_substring = 1
    cdef unsigned int pointer = 0

    # While we haven't decoded the full string we continue
    while prefix_len + len_substring <= stringlength:

        # Given a prefix length, find the largest substring
        if (
            bytestring[pointer + len_substring - 1]
            == bytestring[prefix_len + len_substring - 1]
        ):
            len_substring += 1  # increase the length of the substring
        else:

            max_len_substring = max(len_substring, max_len_substring)
            pointer += 1

            # all the pointers have been investigated, we pick the largest for the jump
            if pointer == prefix_len:
                # Increase the complexity
                complexity += 1

                # Increase the prefix length by the maximum substring size found so far
                prefix_len += max_len_substring

                # Reset the variables
                pointer = 0
                max_len_substring = 1

            # reset the length of the substring
            len_substring = 1

    # Check final repetition if we were in the middle of a substring
    if len_substring != 1:
        complexity += 1

    return complexity

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