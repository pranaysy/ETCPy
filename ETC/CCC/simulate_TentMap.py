#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""

# Import calls
import numpy as np
from numba import vectorize, float64, njit

# Compute single step of iteration through skew-tent map
@vectorize([float64(float64, float64)])
def _skewtent_onestep(value, threshold):
    """
    Computes a single step of iteration through the skew-tent map given an
    input (previous) value and a threshold. Returns the next value as output.
    This function is called by _iterate_skewtent for iterating repeatedly.

    Parameters
    ----------
    value : scalar, float64
        Input value to the skew-tent map.
    threshold : scalar, float64
        Threshold value of the skew-tent map.

    Returns
    -------
    Output value as float64 from the skew-tent map.
    Computed conditionally as follows:
        If value < threshold, then output is value / threshold
        Else, output is (1 - value)/(1 - threshold)

    """
    if value < threshold:
        return value / threshold
    return (1 - value) / (1 - threshold)


# Multiple iterations along skew-tent map
@njit
def _iterate_skewtent(threshold, traj_vec, coupling):
    """
    Computes multiple steps of iteration through the skew-tent map given a
    starting condition, as the first element of an array full of zeros, and
    a threshold for the skew-tent map. This function calls _skewtent_onestep
    for running a single step, and is itself called by _compute_trajectory,
    which initializes the trajectory array.

    Parameters
    ----------
    threshold : vector of size 2, float64
        Threshold value of the skew-tent map.
    traj_vec : array, 2D, float64
        Pre-allocated array of zeroes with the 1st element containing a
        value corresponding to initial condition of the skew-tent map

    Returns
    -------
    traj_vec : array, 2D, float64
        Array populated with values corresponding to the trajectory taken by
        recursive iteration through a skew-tent map. Length of this trajectory
        is inferred from the array shape.

    """
    # Iteration using for-loop over indices
    for idx in range(1, max(traj_vec.shape)):

        # Execute single step of iteration using previous value and threshold
        traj_vec[0, idx] = _skewtent_onestep(traj_vec[0, idx - 1], threshold[0])

        # Linearly dependent tent map
        buffer = _skewtent_onestep(traj_vec[1, idx - 1], threshold[1])
        traj_vec[1, idx] = coupling[0] * traj_vec[0, idx] + (1 - coupling[0]) * buffer

        # Nonlinearly dependent tent map
        buffer = (
            coupling[1] * traj_vec[0, idx - 1]
            + (1 - coupling[1]) * traj_vec[2, idx - 1]
        )
        traj_vec[2, idx] = _skewtent_onestep(buffer, threshold[2])

    # Return populated array
    return traj_vec


# Compute trajectory given initial conditions, threshold and size
@njit
def _compute_trajectory(init_cond, threshold, length, coupling):
    """
    Computes the trajectory along a skew-tent map with given threshold and an
    initial condition for a given distance. Doesn't validate input. This is
    called by compute_trajectory after checking inputs.

    Parameters
    ----------
    init_cond : vector of size 3, float64
        Initial value for iterating through the skew-tent map.
    threshold : vector of size 3, float64
        Threshold value of the skew-tent map.
    length : scalar, integer
        Size of the trajectory to compute through iteration.

    Returns
    -------
    array, 2D, float64
        Array of demanded size filled with values corresponding to the
        trajectory.

    """
    # Pre-allocate array for trajectory with known size
    traj_vec = np.zeros((3, length), dtype=np.float64)

    # Assign initial condition to first elements of Y, Xlin, Xnonlin
    traj_vec[:, 0] = init_cond

    # Run iterations and return populated array
    return _iterate_skewtent(threshold, traj_vec, coupling)


# Warmup for Numba cache initialization
def warmup():
    """
    Runs all the Numba-optimized functions to initialize Numba's JIT.
    Returns nothing and only prints to stdout.

    Returns
    -------
    None.

    """
    initials = np.array([0.1] * 3)
    threshs = np.array([0.2] * 3)
    couplings = np.array([0] * 2)
    expected = np.array([0.625] * 3)
    # Test for a known value
    if (_compute_trajectory(initials, threshs, 3, couplings)[:, -1] == expected).all():
        print("> Numba JIT warmup successful for chaotic_sampler ...")
    else:
        print("> Numba JIT warmup failed for chaotic_sampler ...")


def compute_trajectory(init_cond, threshold, length, burn, coupling):
    """
    Computes the trajectory along a skew-tent map with given threshold and an
    initial condition for a given distance. Wrapper around _compute_trajectory
    and checks inputs for sanity

    Parameters
    ----------
    init_cond : vector of size 3, float64
        Initial value for iterating through the skew-tent map.
            range: 0 < init_cond < 1
    threshold : vector of size 3, float64
        Threshold value of the skew-tent map.
            range: 0 < threshold < 1
    length : scalar, integer
        Size of the trajectory to compute through iteration.
            range: 10^2 < length < 10^7

    Returns
    -------
    array, 2D, float64
        Array of demanded size filled with values corresponding to the
        trajectory.

    """
    # Return trajectory if inputs are valid

    return _compute_trajectory(init_cond, threshold, length + burn, coupling)[:, burn:]


def coupled_TM(threshold, length, burn, coupling, seed):

    np.random.seed(seed)

    # Initialize starting points
    init_cond = np.random.uniform(size=(3))

    # Initialize thresholds
    thresholds = np.array([threshold] * 3)

    # Initialize couplings
    couplings = np.array([coupling] * 2)

    trajectories = compute_trajectory(init_cond, thresholds, length, burn, couplings)

    return {
        "independent": trajectories[0, :],
        "dependent_linear": trajectories[1, :],
        "dependent_nonlinear": trajectories[2, :],
    }
