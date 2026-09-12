"""Scrambled-feature control for the functional-genomics AUC increment.

Tests whether the genetic to genetic+FG AUC increment reflects functional-genomics signal
or added model capacity, by row-permuting the 19 functional-genomics features and
rerunning the increment. See REPRODUCIBILITY.md, Other significance tests, for the
rationale and cost details.

Reuses the CV machinery from attribution_decomposition.py.

Inputs:
  data/perturbseq/pu/pu_model_matrix.parquet

Outputs:
  figure_data/scrambled_feature_control.json
  outputs/permutation/scrambled_feature_control.checkpoint.json (resume state, gitignored)

Run: python scrambled_feature_control.py [--n-scramble 99]
"""
import os, sys, json, argparse
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _label_permutation as lp
import attribution_decomposition as ad

OUT = os.path.join(lp.FIGDATA, "scrambled_feature_control.json")
CKPT = os.path.join(lp.CKPT_DIR, "scrambled_feature_control.checkpoint.json")
FG = ad.PERTURBATIONAL + ad.OBSERV


def main():
    ap = argparse.ArgumentParser()
    # Default n_scramble=99 matches the committed output file; a different default would not match the published figure.
    # 99 puts the rank-based floor at (0+1)/100 = 0.010; draws are seeded RandomState(700+i) and
    # checkpointed individually, so raising this later costs only the additional draws.
    ap.add_argument("--n-scramble", type=int, default=99)
    args = ap.parse_args()

    os.makedirs(lp.CKPT_DIR, exist_ok=True)
    m = pd.read_parquet(ad.MATRIX)
    role = m.pu_role.values
    is_U = role == "unlabeled"
    Pidx, Uidx = np.where(role == "P")[0], np.where(is_U)[0]
    print(f"{len(FG)} functional-genomics features scrambled; "
          f"{len(ad.GENETIC)} genetic features kept intact", file=sys.stderr)

    def incr(mat, seed):
        g = ad.cv_auc(mat, ad.GENETIC, Pidx, Uidx, is_U, seed)
        f = ad.cv_auc(mat, ad.GENETIC + FG, Pidx, Uidx, is_U, seed)
        return np.mean(f) - np.mean(g)

    res = json.load(open(CKPT)) if os.path.exists(CKPT) else {"real": None, "scrambled": {}}

    # Checkpoint state reflects the feature set, depth and N_SEED at run start; delete outputs/permutation/ before rerunning after a config change.
    if res["real"] is None:
        res["real"] = float(np.mean([incr(m, s) for s in range(ad.N_SEED)]))
        json.dump(res, open(CKPT, "w"), indent=1)
        print(f"real increment = {res['real']:+.4f}", file=sys.stderr)

    for i in range(args.n_scramble):
        if str(i) in res["scrambled"]:
            continue
        rng = np.random.RandomState(700 + i)
        mm = m.copy()
        for c in FG:
            mm[c] = mm[c].values[rng.permutation(len(mm))]
        res["scrambled"][str(i)] = float(np.mean([incr(mm, s) for s in range(ad.N_SEED)]))
        json.dump(res, open(CKPT, "w"), indent=1)     # checkpoint each scramble
        print(f"  scramble {i}: {res['scrambled'][str(i)]:+.4f}  "
              f"({len(res['scrambled'])}/{args.n_scramble})", file=sys.stderr)

    # The output JSON is written once at the end so an interrupted run cannot leave a partially populated file; partial state stays in the checkpoint.
    if len(res["scrambled"]) < args.n_scramble:
        raise SystemExit(f"only {len(res['scrambled'])}/{args.n_scramble} scrambles done; "
                         f"rerun to resume from {CKPT}")

    scr = np.array([res["scrambled"][str(i)] for i in range(args.n_scramble)])
    real = res["real"]
    # No z key: with only dozens of draws the permutation SD is too noisy, so emp_p is the cited statistic; n_ge_real, its numerator, is written directly rather than recounted.
    res["summary"] = dict(
        real_increment=real, n_scramble=len(scr),
        scrambled_mean=float(scr.mean()), scrambled_sd=float(scr.std()),
        scrambled_max=float(scr.max()),
        n_ge_real=int(np.sum(scr >= real)),
        emp_p=float((np.sum(scr >= real) + 1) / (len(scr) + 1)),
        note="rank-based emp_p is the statistic to quote; "
             "do NOT compute a z-score from a permutation SD")
    # indent=1 matches the committed output file's byte layout.
    json.dump(res, open(OUT, "w"), indent=1)
    print(json.dumps(res["summary"], indent=1))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
