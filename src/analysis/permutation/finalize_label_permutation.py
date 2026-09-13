"""Recompute label_permutation_pvalues.csv from an existing null table.
See REPRODUCIBILITY.md.

Usage: python finalize_label_permutation.py
"""
import os, sys, json
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _label_permutation as lp

CKPT = os.path.join(lp.CKPT_DIR, "label_permutation_null.checkpoint.json")


def load_null():
    """Prefers a checkpoint from a live run when present; falls back to the committed null table otherwise."""
    if os.path.exists(CKPT):
        done = {int(k): v for k, v in json.load(open(CKPT)).items()}
        null = pd.DataFrame([{"perm": i, **done[i]} for i in sorted(done)])
        return null, f"checkpoint ({len(null)} perms)"
    if os.path.exists(lp.NULL_CSV):
        # No check.names equivalent: rung columns contain '+', so read them verbatim to avoid mismatches.
        return pd.read_csv(lp.NULL_CSV), f"committed {os.path.basename(lp.NULL_CSV)}"
    raise SystemExit(
        f"no null table to finalize from.\n  looked for: {CKPT}\n  and:        {lp.NULL_CSV}\n"
        "Run `make permutation-null` to generate one (~2 h on 11 threads).")


def main():
    null, src = load_null()
    obs = lp.observed_ladder()

    prev = pd.read_csv(lp.PVAL_CSV) if os.path.exists(lp.PVAL_CSV) else None
    out = lp.pvalue_table(null, obs)
    out.to_csv(lp.PVAL_CSV, index=False)

    print(f"null source     : {src}")
    print(f"observed ladder : {lp.LADDER_CSV}")
    for k in lp.RUNGS:
        print(f"  {k:24s} {obs[k]!r}")

    if prev is not None:
        merged = prev.merge(out, on=["metric", "kind"], suffixes=("_old", "_new"))
        changed = merged[(merged.observed_old != merged.observed_new) |
                         (merged.emp_p_old != merged.emp_p_new)]
        if changed.empty:
            print("\nno change: p-values were already current against this ladder")
        else:
            print(f"\n{len(changed)} metric(s) changed:")
            for _, r in changed.iterrows():
                print(f"  {r.metric} ({r.kind})")
                print(f"    observed {r.observed_old!r} -> {r.observed_new!r}")
                if r.emp_p_old != r.emp_p_new:
                    print(f"    emp_p    {r.emp_p_old!r} -> {r.emp_p_new!r}")
    print(f"\nwrote {lp.PVAL_CSV}")


if __name__ == "__main__":
    main()
