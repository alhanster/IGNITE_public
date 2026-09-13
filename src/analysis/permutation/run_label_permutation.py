"""Label-permutation null for the evidence AUC ladder (GPS-style, Duffy/Do 2024).

Opt-in: run via `make permutation-null`, which checkpoints every 25 permutations and resumes;
see REPRODUCIBILITY.md, permutation-null / permutation-null-finalize.

Inputs (read-only):
  data/perturbseq/pu/pu_model_matrix.parquet   labels and features
  figure_data/panelA_auc_ladder.csv            observed ladder

Outputs:
  figure_data/label_permutation_null.csv       per-permutation null AUC per rung and increment
  figure_data/label_permutation_pvalues.csv    observed value, null mean/sd, empirical p
  outputs/permutation/label_permutation_null.checkpoint.json   resume state (gitignored)

Run: python run_label_permutation.py [--n-perm 1000] [--jobs 11]
"""
import os, sys, json, argparse, time
# Thread env vars set before numpy/xgboost import; later imports pull in pandas and xgboost.
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
from concurrent.futures import ThreadPoolExecutor
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "model"))
import _label_permutation as lp
import attribution_decomposition as ad   # unmodified CV machinery; also runs check_pins()
from xgboost import XGBClassifier

# Pins each base learner to a single XGBoost thread; see REPRODUCIBILITY.md, permutation-null / permutation-null-finalize.
_MAX_DEPTH = int(os.environ.get("MODEL_MAX_DEPTH", "2"))


def _single_thread_learner(seed):
    return XGBClassifier(n_estimators=120, max_depth=_MAX_DEPTH, learning_rate=0.1,
                         subsample=0.8, colsample_bytree=0.8, eval_metric="logloss",
                         n_jobs=1, verbosity=0, random_state=seed)


ad.base_learner = _single_thread_learner

CKPT = os.path.join(lp.CKPT_DIR, "label_permutation_null.checkpoint.json")
CKPT_EVERY = 25
RUNG_COLS = [("genetic", ad.GENETIC),
             ("genetic+observational", ad.GENETIC + ad.OBSERV),
             ("genetic+perturbational", ad.GENETIC + ad.PERTURBATIONAL),
             ("full", ad.GENETIC + ad.PERTURBATIONAL + ad.OBSERV)]


def one_perm(i, m, pool, n_P):
    """One permutation: reassigns the positive label to n_P random genes drawn from the scorable
    pool, then runs the ladder.

    m and pool are shared read-only across threads; each call builds its own label masks
    without mutating them. Trial-heldout genes are excluded from pool, so a permutation never
    assigns the positive label to them.
    """
    rng = np.random.RandomState(10_000 + i)          # deterministic per-perm seed
    fakeP = rng.choice(pool, size=n_P, replace=False)
    is_fakeP = np.zeros(len(m), bool); is_fakeP[fakeP] = True
    is_fakeU = np.zeros(len(m), bool); is_fakeU[pool] = True; is_fakeU[fakeP] = False
    Pi, Ui = np.where(is_fakeP)[0], np.where(is_fakeU)[0]
    out = {}
    for name, cols in RUNG_COLS:
        aucs = []
        for s in range(ad.N_SEED):
            aucs += ad.cv_auc(m, cols, Pi, Ui, is_fakeU, seed=s)
        out[name] = float(np.mean(aucs))
    return i, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--jobs", type=int, default=max(1, os.cpu_count() - 1))
    args = ap.parse_args()

    os.makedirs(lp.CKPT_DIR, exist_ok=True)
    m = pd.read_parquet(ad.MATRIX)
    role = m.pu_role.values
    is_P, is_U = role == "P", role == "unlabeled"
    n_P = int(is_P.sum())
    pool = np.where(is_P | is_U)[0]                   # trial_heldout excluded
    obs = lp.observed_ladder()
    print(f"genes={len(m)} P={n_P} scorable_pool={len(pool)} "
          f"heldout={int((role == 'trial_heldout').sum())}", file=sys.stderr)
    print(f"observed ladder: {obs}", file=sys.stderr)

    done = {}
    if os.path.exists(CKPT):
        done = {int(k): v for k, v in json.load(open(CKPT)).items()}
        print(f"resuming: {len(done)} perms already done", file=sys.stderr)

    todo = [i for i in range(args.n_perm) if i not in done]
    t0, n_start = time.time(), len(done)
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        for b0 in range(0, len(todo), CKPT_EVERY):
            batch = todo[b0:b0 + CKPT_EVERY]
            for f in [ex.submit(one_perm, i, m, pool, n_P) for i in batch]:
                i, r = f.result()
                done[i] = r
            json.dump({str(k): v for k, v in done.items()}, open(CKPT, "w"))
            el, n_new = time.time() - t0, len(done) - n_start
            print(f"  {len(done)}/{args.n_perm} perms  ({el:.0f}s, "
                  f"{el / max(1, n_new):.1f}s/perm)", file=sys.stderr)

    null = pd.DataFrame([{"perm": i, **done[i]} for i in sorted(done)])
    null.to_csv(lp.NULL_CSV, index=False)
    lp.pvalue_table(null, obs).to_csv(lp.PVAL_CSV, index=False)
    print(f"wrote {lp.NULL_CSV} and {lp.PVAL_CSV} ({len(null)} perms)", file=sys.stderr)


if __name__ == "__main__":
    main()
