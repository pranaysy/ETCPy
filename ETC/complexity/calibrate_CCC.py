#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""

import ETC
from ETC import compute_1D, compute_2D, partition
from itertools import product
from collections import defaultdict
from functools import partial
from multiprocessing import Pool
from time import perf_counter
from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns
sns.set()

get1D = partial(compute_1D, order=2, verbose=False, truncate=True)
get2D = partial(compute_2D, order=2, verbose=False, truncate=True)

def test(seq_x, seq_y, past_win_size, delta, step_size, partitions=False):
    if partitions:
        seq_x = partition(seq_x, partitions)
        seq_y = partition(seq_y, partitions)

    out = defaultdict(list)

    length = len(seq_x)
    total_win_size = past_win_size + delta

    for n, k in enumerate(range(0, length - total_win_size, step_size)):

        out["window"].append(n+1)
        out["begin"].append(k)
        out["past_win_size"].append(past_win_size)
        out["end_past"].append(k + past_win_size)
        out["delta"].append(delta)
        out["total_win_size"].append(total_win_size)
        out["end_total"].append(k + total_win_size)
        out["step_size"].append(step_size)

        ## CC 1D for X --------------------------------------------------------
        # ETC 1D for past values of X
        ETC1D_X_ini = get1D(seq_x[k : k + past_win_size])["ETC1D"]
        out["ETC_1D_X_past_raw"].append(ETC1D_X_ini)

        ETC1D_X_ini /= past_win_size - 1
        out["ETC_1D_X_past_norm"].append(ETC1D_X_ini)

        # ETC 1D for past and current+past=total values of X
        ETC1D_X_fin = get1D(seq_x[k : k + total_win_size])["ETC1D"]
        out["ETC_1D_X_total_raw"].append(ETC1D_X_fin)

        ETC1D_X_fin /= total_win_size - 1
        out["ETC_1D_X_total_norm"].append(ETC1D_X_fin)

        # CC 1D for past and total values of X
        CC1D_X_past = ETC1D_X_fin - ETC1D_X_ini
        out["CC_1D_X"].append(CC1D_X_past)
        ## --------------------------------------------------------------------

        ## CC 1D for Y --------------------------------------------------------
        # ETC 1D for past values of Y
        ETC1D_Y_ini = get1D(seq_y[k : k + past_win_size])["ETC1D"]
        out["ETC_1D_Y_past_raw"].append(ETC1D_Y_ini)

        ETC1D_Y_ini /= past_win_size - 1
        out["ETC_1D_Y_past_norm"].append(ETC1D_Y_ini)

        # ETC 1D for past and current+past=total values of X
        ETC1D_Y_fin = get1D(seq_y[k : k + total_win_size])["ETC1D"]
        out["ETC_1D_Y_total_raw"].append(ETC1D_Y_fin)

        ETC1D_Y_fin /= total_win_size - 1
        out["ETC_1D_Y_total_norm"].append(ETC1D_Y_fin)

        # CC 1D for past and total values of X
        CC1D_Y_past = ETC1D_Y_fin - ETC1D_Y_ini
        out["CC_1D_Y"].append(CC1D_Y_past)
        ## --------------------------------------------------------------------

        # ETC 2D for past values of X and Y -----------------------------------
        ETC2D_ini = get2D(seq_x[k : k + past_win_size], seq_y[k : k + past_win_size])[
            "ETC2D"
        ]
        out["ETC_2D_X_past_Y_past_raw"].append(ETC2D_ini)
        out["ETC_2D_Y_past_X_past_raw"].append(ETC2D_ini)

        ETC2D_ini /= past_win_size - 1
        out["ETC_2D_X_past_Y_past_norm"].append(ETC2D_ini)
        out["ETC_2D_Y_past_X_past_norm"].append(ETC2D_ini)
        ## --------------------------------------------------------------------

        # ETC 2D for current+past=total values of X and past values of Y plus current values of X
        ETC2D_X_fin = get2D(
            seq_x[k : k + total_win_size],
            seq_y[k : k + past_win_size]
            + seq_x[k + past_win_size : k + total_win_size],
        )["ETC2D"]
        out["ETC_2D_X_total_Y_past_raw"].append(ETC2D_X_fin)

        ETC2D_X_fin /= total_win_size - 1
        out["ETC_2D_X_total_Y_past_norm"].append(ETC2D_X_fin)

        # CC 2D for past and total values of X
        CC2D_X_total_Y_past = ETC2D_X_fin - ETC2D_ini
        out["CC_2D_X_by_Y_past"].append(CC2D_X_total_Y_past)
        ## --------------------------------------------------------------------

        # ETC 2D for current+past=total values of Y and past values of X plus current values of Y
        ETC2D_Y_fin = get2D(
            seq_y[k : k + total_win_size],
            seq_x[k : k + past_win_size]
            + seq_y[k + past_win_size : k + total_win_size],
        )["ETC2D"]
        out["ETC_2D_Y_total_X_past_raw"].append(ETC2D_Y_fin)

        ETC2D_Y_fin /= total_win_size - 1
        out["ETC_2D_Y_total_X_past_norm"].append(ETC2D_Y_fin)

        # CC 2D for past and total values of X
        CC2D_Y_total_X_past = ETC2D_Y_fin - ETC2D_ini
        out["CC_2D_Y_by_X_past"].append(CC2D_Y_total_X_past)
        ## --------------------------------------------------------------------

    return pd.DataFrame(out)

def test_multiple(seq_x, seq_y):

    # Past window size
    PWS = [100,150,175,200]

    # Current window size
    CWS = [10,15,20,25]

    # Jump step size
    SS = [10,15,20,25,30]

    before = perf_counter()

    results = []

    for past_win_size, delta, step_size in product(PWS, CWS,SS):
        results.append(test(seq_x, seq_y, past_win_size, delta, step_size))

    results = pd.concat(results)

    after = perf_counter()

    return results, after - before

def apply_func(params):
    past_win_size, delta, step_size = params
    return test(seq_x, seq_y, past_win_size, delta, step_size)

def test_multiple_parallel(seq_x, seq_y):

    # Past window size
    PWS = range(10,201,10)

    # Current window size
    CWS = [15]

    # Jump step size
    SS = range(10,201,50)

    before = perf_counter()
    # Initialize pool of parallel workers
    pool = Pool()

    # Map-execute function across files
    results = pool.map(apply_func, product(PWS, CWS,SS))

    # Graceful exit
    pool.close()
    pool.join()

    results = pd.concat(results)

    after = perf_counter()

    return results, after - before

a2, timings = test_multiple_parallel(seq_x, seq_y)


# %%
fig, ax = plt.subplots(1,1)
sns.lineplot(data=a2, x='past_win_size', y='ETC_1D_X_past_norm', ax=ax)
sns.lineplot(data=a2, x='past_win_size', y='ETC_1D_X_total_norm', ax=ax)
sns.lineplot(data=a2, x='past_win_size', y='ETC_2D_X_past_Y_past_norm', ax=ax)
sns.lineplot(data=a2, x='past_win_size', y='ETC_2D_X_total_Y_past_norm', ax=ax)