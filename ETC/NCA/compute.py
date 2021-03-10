#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module has functions for multicore computation of NCA from a 2D matrix

@author: Pranay S. Yadav
"""
import pandas as pd
from ETC.NCA import parallelize_jl as NCAP


def compute_CCC(matrix, CCC_params):
    """
    Compute causal complexity estimates for all pairs of rows of input matrix

    Estimates are derived as per 4 models: ETCP, ETCE, LZP and CCC

    Parameters
    ----------
    matrix : np.ndarray, 2d, uint32
        MxN matrix, C(M,2) rowpairs, with each row of length N.
    CCC_params : dict
        CCC parameters with the following names for keys:
            "LEN_past", "ADD_meas", "STEP_size"

    Returns
    -------
    pd.DataFrame
        DataFrame containing causal estimates from all 4 models.

    """
    # Create a generator that produces rowpairs one at a time
    rowpairs = NCAP.get_rowpairs(matrix)

    # Compute causal estimates in parallel across rowpairs
    estimates = NCAP.parallelized_CCC(rowpairs, CCC_params)

    # Convert estimates to a DataFrame and return
    return pd.DataFrame(estimates)


def compute_CCM(matrix, kernel="LZ"):
    """
    Compute causal complexity estimates for all pairs of rows of input matrix

    Estimates are derived as per 4 models: ETCP, ETCE, LZP and CCC

    Parameters
    ----------
    matrix : np.ndarray, 2d, uint32
        MxN matrix, C(M,2) rowpairs, with each row of length N.
    CCC_params : dict
        CCC parameters with the following names for keys:
            "LEN_past", "ADD_meas", "STEP_size"

    Returns
    -------
    pd.DataFrame
        DataFrame containing causal estimates from all 4 models.

    """
    # Create a generator that produces rowpairs one at a time
    rowpairs = NCAP.get_rowpairs(matrix)

    # Compute causal estimates in parallel across rowpairs
    estimates = NCAP.parallelized_CCM(rowpairs, kernel)

    # Convert estimates to a DataFrame and return
    return pd.DataFrame(estimates)


def get_causal(df):
    """
    Extract causal strengths from all estimates as a long-form (tidy) DataFrame

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing raw causal estimates from all 4 models.
        As returned by get_estimates()

    Returns
    -------
    pd.DataFrame
        DataFrame containing only causal strengths from the 4 models in both directions.

    """
    # Fix pair identifiers
    identifiers = ["index_pair", "index_x", "index_y"]

    # Initialize aggregator of melted / tidied dataframes
    dfs = list()

    # Iterate over each model and tidy up
    for model in ["ETCP", "ETCE", "LZP", "CCC"]:

        if df.filter(like=model).shape[-1] != 0:
            dat = df.melt(
                id_vars=identifiers,
                value_vars=[f"{model}_x_to_y", f"{model}_y_to_x"],
                var_name="direction",
                value_name=model,
            )

            # Strip column names to make it generic -> easier to concatenate
            dat["direction"] = dat["direction"].str.replace(f"{model}_", "")

            # Set identifiers as indices and append to aggregator
            dfs.append(dat.set_index(identifiers + ["direction"]))

    # Concatenate all tidied dataframes and return
    return pd.concat(dfs, axis=1).reset_index()


def get_NCA(df, k=0.1):
    """
    Compute NCA from top k causal strengths (top k pairs)

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing only causal strengths from the 4 models in both directions
        in tidy/long-form, as returned by get_causal()

    k : float, optional, 0 < k < 1
        Top proportion of causal strengths. The default is 0.1.

    Returns
    -------
    pd.DataFrame
        DataFrame containing NCA estimates with summary statistics.

    """
    # Convert k to numeric based on available causal pairs
    k_int = round(k * len(df))

    # Initialize aggregator of estimates from each model
    agg = list()

    # Iterate over each model, compute summary statistics and aggregate
    for model in ["ETCP", "ETCE", "LZP", "CCC"]:

        if df.filter(like=model).shape[-1] != 0:
            dat = df[model].nlargest(k_int)
            agg.append(
                {
                    "model": model,
                    "mean": dat.mean(),
                    "median": dat.median(),
                    "max": dat.max(),
                    "min": dat.min(),
                    "mad": dat.mad(),
                    "std": dat.std(),
                }
            )

    # Combine all estimates into a DataFrame and return
    return pd.DataFrame(agg).set_index("model")


# Function for earlier versions, preserved for later
# def get_complexity(df):
#     """
#     Extract CCM estimates of complexity - ETC and LZ of each row element

#     Parameters
#     ----------
#     df : pd.DataFrame
#         DataFrame containing raw causal estimates from all 4 models.
#         As returned by get_estimates()

#     Returns
#     -------
#     pd.DataFrame
#         DataFrame containing unique ETC and LZ estimates for each row.

#     """
#     # Store estimates for one row of the pair
#     dfx = df[["index_x", "ETC_x", "LZ_x"]]
#     dfx = dfx.rename(columns=lambda m: m.replace("_x", ""))

#     # Store estimates for the other row of the pair
#     dfy = df[["index_y", "ETC_y", "LZ_y"]]
#     dfy = dfy.rename(columns=lambda m: m.replace("_y", ""))

#     # Combine the two, drop duplicates and return
#     return pd.concat([dfx, dfy]).drop_duplicates().set_index("index")
