#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""
import numpy as np
import pandas as pd
from ETC.NCA import parallelize as NCAP


def _get_estimates(matrix, CCC_params):

    rowpairs = NCAP.get_rowpairs(matrix)

    estimates = NCAP.parallelized(rowpairs, CCC_params)
    df = pd.DataFrame(estimates)

    return df


def get_causal(df):

    identifiers = ["index_pair", "index_x", "index_y"]

    dfs = list()
    for model in ["ETCP", "ETCE", "LZP", "CCC"]:
        dat = df.melt(
            id_vars=identifiers,
            value_vars=[f"{model}_x_to_y", f"{model}_y_to_x"],
            var_name="direction",
            value_name=model,
        )
        dat["direction"] = dat["direction"].str.replace(f"{model}_", "")
        dfs.append(dat.set_index(identifiers + ["direction"]))

    return pd.concat(dfs, axis=1).reset_index()


def get_NCA(df, k=0.1):

    k_int = round(k * len(df))

    agg = list()
    for model in ["ETCP", "ETCE", "LZP", "CCC"]:
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

    return pd.DataFrame(agg).set_index("model")


def get_complexity(df):

    dfx = df[["index_x", "ETC_x", "LZ_x"]]
    dfx = dfx.rename(columns=lambda m: m.replace("_x", ""))

    dfy = df[["index_y", "ETC_y", "LZ_y"]]
    dfy = dfy.rename(columns=lambda m: m.replace("_y", ""))

    return pd.concat([dfx, dfy]).drop_duplicates().set_index("index")


# %%
matr = np.random.randint(1, 3, [10, 5000], dtype="uint32")

q = {"LEN_past": 250, "ADD_meas": 500, "STEP_size": 250}
out = _get_estimates(matr, q)

# %%
# a = NCAP._kernel((1, (0, 2, matr[0, :], matr[1, :])), q)
