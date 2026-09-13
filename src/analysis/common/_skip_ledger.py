#!/usr/bin/env python3
import os

LEDGER = os.path.join(
    os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..")),
    "outputs", "_skipped.tsv")


def record_skip(script, reason, outputs):
    """Append one skip record.

    script  : basename of the step that skipped
    reason  : short phrase naming the absent input
    outputs : iterable of figure_data paths left at their committed values
    """
    os.makedirs(os.path.dirname(LEDGER), exist_ok=True)
    new = not os.path.exists(LEDGER)
    with open(LEDGER, "a", encoding="utf-8") as fh:
        if new:
            fh.write("script\treason\toutputs_left_committed\n")
        fh.write("%s\t%s\t%s\n" % (script, reason, ";".join(outputs)))


def read_skips():
    """Return the recorded skips as a list of dicts, empty when nothing skipped."""
    if not os.path.exists(LEDGER):
        return []
    rows = []
    with open(LEDGER, encoding="utf-8") as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip()]
    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) >= 3:
            rows.append({"script": parts[0], "reason": parts[1],
                         "outputs": [p for p in parts[2].split(";") if p]})
    return rows


if __name__ == "__main__":
    skips = read_skips()
    if not skips:
        # Deliberately not "every figure_data table": the 11 opt-in tables (discordance_*,
        # label_permutation_*, vignette_panelC_gwas_*) have no producer in `make tables` at
        # all, so they are not skips this ledger can see, and they were not re-derived either.
        print("No steps skipped: every table `make tables` produces was re-derived on this run.")
        print("  The 11 opt-in tables have no producer in this stage and are unchanged;")
        print("  see REPRODUCIBILITY.md, Reproduction boundary.")
        raise SystemExit(0)
    n_out = sum(len(s["outputs"]) for s in skips)
    print("")
    print("=" * 72)
    print("%d step(s) skipped; %d figure_data table(s) were NOT re-derived on this run"
          % (len(skips), n_out))
    print("=" * 72)
    for s in skips:
        print("  %s" % s["script"])
        print("      reason: %s" % s["reason"])
        for o in s["outputs"]:
            print("      committed value used: %s" % o)
    print("")
    print("  These tables are the committed values, not products of this run, so")
    print("  `make verify-tables` compares them against themselves and passes trivially.")
    print("  See REPRODUCIBILITY.md, Reproduction boundary.")
    print("")
