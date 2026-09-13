"""See REPRODUCIBILITY.md."""
from scipy.stats import false_discovery_control


def holm(ps):
    """Holm-Bonferroni step-down adjusted p-values."""
    n = len(ps)
    order = sorted(range(n), key=lambda i: ps[i])
    adj = [0.0] * n
    run = 0.0
    for rank, i in enumerate(order):
        run = max(run, min(1.0, float(ps[i]) * (n - rank)))
        adj[i] = run
    return adj


def bh(ps):
    """Benjamini-Hochberg adjusted p-values, in the same order as the input; callers must preserve correspondence to labels themselves."""
    return [float(v) for v in false_discovery_control(list(ps), method="bh")]
