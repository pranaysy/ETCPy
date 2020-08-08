#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""

import numpy as np


def coupled_AR(length=1000, a=0.9, b=0.8, c=0.8, e=0.01, burn=100, seed=1):
    """
    Generate discrete-time coupled AR processes with known parameters.

    Dependent process defined as:
        x[n] = a * x[n - 1] + b * y[n - 1] + e * noise_x[n]

    Independent process defined as:
        y[n] = c * y[n - 1] + e * noise_y[n]

    Parameters
    ----------
    length : int, optional
        Legnth of samples drawn from the process. The default is 1000.
    a : float, optional
        Coefficient for dependent process, capturing dependency on its own past.
        The default is 0.9.
    b : float, optional
        Coefficient for dependent process, capturing dependency on the independent
        process - causal interaction from independent to dependent. The default is 0.8.
    c : float, optional
        Coefficient for independent process, capturing dependency on its own past.
        The default is 0.8.
    e : float, optional
        Coefficient for uniform random noise mixture. The default is 0.01
    burn : int, optional
        Number of initial samples to burn. The default is 100.
    seed: int, optional
        Seed value for initialization of random number generator. The default is 1

    Returns
    -------
    dict
        Two key-value pairs:
            "dependent": Samples of the dependent process.
            "independent": Samples of the independent process.

    """
    # Anchor seed for reproducibility
    np.random.seed(seed)

    # AR processes: initialize
    x = np.zeros(length, dtype="float64")
    y = np.zeros(length, dtype="float64")

    # Generate noise vector of appropriate length
    noise_x = e * np.random.normal(0, 1, length + burn)
    noise_y = e * np.random.normal(0, 1, length + burn)

    # Initialize starting points
    x[0] = np.random.uniform()
    y[0] = np.random.uniform()

    # Burn initial samples
    if burn:
        for n in range(burn):
            x[0] = a * x[0] + b * y[0] + noise_x[n]
            y[0] = c * y[0] + noise_y[n]

    # Store further samples
    for n in range(1, length):
        x[n] = a * x[n - 1] + b * y[n - 1] + noise_x[n]
        y[n] = c * y[n - 1] + noise_y[n]

    return {"dependent": x, "independent": y}
