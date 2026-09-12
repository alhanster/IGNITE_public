#!/usr/bin/env python
"""See REPRODUCIBILITY.md, Other significance tests."""
import numpy as np
import pandas as pd


def hits_at_k(rank_pctile, is_emergent, ks):
    """rank_pctile over the pool (higher is better); is_emergent is boolean and index-aligned, used for top-k selection."""
    order = np.argsort(-rank_pctile)
    out = {}
    for k in ks:
        topk = order[:k]
        out[k] = int(is_emergent[topk].sum())
    return out


def recovery_with_ci(scores_df, score_col, ks=(50, 100, 200, 500),
                     n_boot=2000, n_perm=10000, seed=0):
    """scores_df requires in_pool, is_emergent, and score_col. Ranks are computed within the pool defined by in_pool, not globally across rows."""
    rng = np.random.RandomState(seed)
    pool = scores_df[scores_df["in_pool"]].copy()
    s = pool[score_col].values.astype(float)
    e = pool["is_emergent"].values.astype(bool)
    # rank within pool; nan scores get the lowest rank
    s_filled = np.where(np.isnan(s), -np.inf, s)
    n_emg = int(e.sum()); n_pool = len(pool)

    obs = hits_at_k(s_filled, e, ks)

    # gene-level bootstrap CI: resample pool rows with replacement
    boot = {k: np.empty(n_boot) for k in ks}
    for b in range(n_boot):
        idx = rng.randint(0, n_pool, n_pool)
        hb = hits_at_k(s_filled[idx], e[idx], ks)
        for k in ks: boot[k][b] = hb[k]

    # permutation null: ranking fixed, emergent labels shuffled over pool
    null = {k: np.empty(n_perm) for k in ks}
    order = np.argsort(-s_filled)
    ranks_of_topk = {k: set(order[:k].tolist()) for k in ks}
    for p in range(n_perm):
        perm_idx = rng.choice(n_pool, n_emg, replace=False)
        pe = np.zeros(n_pool, bool); pe[perm_idx] = True
        for k in ks:
            null[k][p] = pe[order[:k]].sum()

    rows = []
    for k in ks:
        nk = null[k]
        exp = nk.mean()
        p_emp = (np.sum(nk >= obs[k]) + 1) / (n_perm + 1)
        rows.append(dict(
            k=k, hits=obs[k],
            hits_ci95_lo=np.percentile(boot[k], 2.5),
            hits_ci95_hi=np.percentile(boot[k], 97.5),
            expected_null=exp,
            null_ci95_lo=np.percentile(nk, 2.5),
            null_ci95_hi=np.percentile(nk, 97.5),
            enrichment=obs[k]/exp if exp > 0 else np.nan,
            p_emp=p_emp,
            n_emergent=n_emg, n_pool=n_pool,
        ))
    return pd.DataFrame(rows)


def paired_bootstrap_delta(scores_df, col_a, col_b, ks=(50,100,200,500), n_boot=2000, seed=0):
    """Paired gene-level bootstrap of hits_a minus hits_b (e.g. full vs genetics model) with CI and p-value; the same resample is applied to both."""
    rng = np.random.RandomState(seed)
    pool = scores_df[scores_df["in_pool"]].copy()
    e = pool["is_emergent"].values.astype(bool)
    sa = np.where(np.isnan(pool[col_a].values), -np.inf, pool[col_a].values)
    sb = np.where(np.isnan(pool[col_b].values), -np.inf, pool[col_b].values)
    n = len(pool)
    obs = {k: hits_at_k(sa,e,[k])[k]-hits_at_k(sb,e,[k])[k] for k in ks}
    boot = {k: np.empty(n_boot) for k in ks}
    for b in range(n_boot):
        idx = rng.randint(0,n,n)
        for k in ks:
            boot[k][b] = hits_at_k(sa[idx],e[idx],[k])[k]-hits_at_k(sb[idx],e[idx],[k])[k]
    rows=[]
    for k in ks:
        bk=boot[k]
        p_two = 2*min((bk<=0).mean(),(bk>=0).mean())
        rows.append(dict(k=k, delta_full_minus_genetics=obs[k],
                         delta_ci95_lo=np.percentile(bk,2.5),
                         delta_ci95_hi=np.percentile(bk,97.5),
                         p_two_sided=min(p_two,1.0)))
    return pd.DataFrame(rows)
