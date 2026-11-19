"""Complete combinatorial measurement of epistasis in Cas12j engineering.

Reproduces every statistic of the combinatorial-library analysis reported in
"Complete combinatorial mutagenesis resolves antagonistic epistasis and yields
a hypercompact Cas12j genome editor".

Inputs (data/)
    combinatorial_library_N.csv   N-segment library, 2,048 designed genotypes
    combinatorial_library_C.csv   C-segment library, 1,024 designed genotypes
    substitution_assignment.csv   library column -> substitution name

Outputs (results/)
    library_summary.csv                      filtering, replicate diagnostics
    pairwise_interaction_detectability.csv   per-pair interaction values
    interaction_vs_background_load.csv       mean interaction by background load
    local_optima.csv                         optima and basins of attraction
    substitution_background_effects.csv      per-substitution contributions
    strategy_benchmark.csv                   in silico strategy comparison
    model_order_cv.csv                       cross-validated model order
    coefficient_transfer.csv                 low-order vs complete coefficients

Usage:  python epistasis_analysis.py [--data DIR] [--out DIR]

Python 3.13, numpy, pandas, scipy.  All random draws use the seeds below, so
the simulation results are reproducible exactly.
"""

import argparse
import itertools
import os

import numpy as np
import pandas as pd
from numpy.linalg import lstsq
from scipy import stats

# --------------------------------------------------------------------------
# constants fixed by the experiment
# --------------------------------------------------------------------------
READ_COUNT_THRESHOLD = {"N": 888, "C": 1647}

# Standard error of a single genotype score, estimated from replicate
# agreement in the primary analysis and used for every error propagation.
# library_summary.csv reports the replicate diagnostics recomputed here.
SCORE_SE = {"N": 0.092, "C": 0.031}

SEED_STRATEGY = 3          # in silico strategy comparison
SEED_CV = 0                # cross-validation folds
N_STRATEGY_RUNS = 200
N_VALIDATION = 10          # candidates re-measured before selection
CV_FOLDS = 10


# --------------------------------------------------------------------------
# input
# --------------------------------------------------------------------------
def load_libraries(data_dir):
    """Return {segment: dict(X, y, labels, replicates, n_designed)}."""
    asg = pd.read_csv(os.path.join(data_dir, "substitution_assignment.csv"))
    spec = {
        "N": dict(path="combinatorial_library_N.csv", seq="AA seq",
                  count="count_data", score="mean_all",
                  reps=["reporter1_rep1", "reporter1_rep2",
                        "reporter2_rep1", "reporter2_rep2"]),
        "C": dict(path="combinatorial_library_C.csv", seq="AAseq",
                  count="total read count", score="mean",
                  reps=["replication1", "replication2"]),
    }
    out = {}
    for seg, s in spec.items():
        raw = pd.read_csv(os.path.join(data_dir, s["path"]))
        kept = raw[raw[s["count"]] >= READ_COUNT_THRESHOLD[seg]].reset_index(drop=True)
        labels = (asg[asg.library == seg].sort_values("col")["substitution"]
                  .tolist())
        mutant_aa = [lab[-1] for lab in labels]
        X = np.array([[1 if seq[i] == aa else 0 for i, aa in enumerate(mutant_aa)]
                      for seq in kept[s["seq"]]], dtype=int)
        out[seg] = dict(X=X,
                        y=kept[s["score"]].values.astype(float),
                        labels=labels,
                        replicates=kept[s["reps"]].values.astype(float),
                        n_designed=len(raw))
    return out


def library_summary(lib):
    """Filtering and replicate diagnostics, including candidate SE estimators."""
    rows = []
    for seg, d in lib.items():
        R = d["replicates"]
        k = R.shape[1]
        pooled_within = float(np.sqrt(np.mean(np.var(R, axis=1, ddof=1))))
        pairs = [R[:, a] - R[:, b] for a, b in itertools.combinations(range(k), 2)]
        paired = float(np.std(np.concatenate(pairs), ddof=1) / np.sqrt(2))
        rows.append(dict(library=f"{seg}-segment",
                         genotypes_designed=d["n_designed"],
                         genotypes_retained=len(d["y"]),
                         retained_pct=round(100 * len(d["y"]) / d["n_designed"], 1),
                         read_count_threshold=READ_COUNT_THRESHOLD[seg],
                         substitutions=len(d["labels"]),
                         replicate_measurements=k,
                         se_used=SCORE_SE[seg],
                         se_pooled_within_over_sqrt_k=round(pooled_within / np.sqrt(k), 4),
                         se_paired_differences_over_sqrt_k=round(paired / np.sqrt(k), 4),
                         best_score=round(float(d["y"].max()), 4),
                         substitutions_in_best=int(d["X"][int(np.argmax(d["y"]))].sum())))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# model-free pairwise epistasis
# --------------------------------------------------------------------------
def contrasts(X, y):
    """All four-genotype contrasts: {(i, j): [(load, epsilon), ...]}."""
    L = X.shape[1]
    index = {tuple(row): n for n, row in enumerate(map(tuple, X))}
    out = {}
    for i, j in itertools.combinations(range(L), 2):
        others = [k for k in range(L) if k not in (i, j)]
        rec = []
        for bits in itertools.product([0, 1], repeat=len(others)):
            g = [0] * L
            for k, b in zip(others, bits):
                g[k] = b
            quad = []
            for a, b in ((0, 0), (1, 0), (0, 1), (1, 1)):
                h = g.copy()
                h[i], h[j] = a, b
                quad.append(index.get(tuple(h)))
            if any(q is None for q in quad):
                continue
            eps = y[quad[3]] - y[quad[1]] - y[quad[2]] + y[quad[0]]
            rec.append((int(sum(bits)), float(eps)))
        out[(i, j)] = rec
    return out


def pairwise_table(lib):
    """Per-pair interaction values, detection limits and significance."""
    frames = []
    for seg, d in lib.items():
        se = SCORE_SE[seg]
        cycle_se = 2 * se                      # contrast of four measurements
        rows = []
        for (i, j), rec in contrasts(d["X"], d["y"]).items():
            loads = np.array([r[0] for r in rec])
            eps = np.array([r[1] for r in rec])
            m = len(eps)
            mean = float(eps.mean())
            parental = float(eps[loads == 0][0]) if (loads == 0).any() else np.nan
            sd_bg = float(eps.std(ddof=1))
            # does the spread across backgrounds exceed measurement noise?
            chi2 = (m - 1) * sd_bg ** 2 / cycle_se ** 2
            p_excess = float(stats.chi2.sf(chi2, m - 1))
            rows.append(dict(pair=f"{i + 1}-{j + 1}", n_bg=m,
                             eps_mean=round(mean, 5),
                             eps_wt=round(parental, 5),
                             sig_complete=bool(abs(mean) > 2 * cycle_se / np.sqrt(m)),
                             sig_single_bg=bool(abs(parental) > 2 * cycle_se),
                             sd_across_bg=round(sd_bg, 5),
                             excess_sd=round(float(max(0.0, sd_bg ** 2 - cycle_se ** 2) ** 0.5), 5),
                             p_excess=p_excess))
        t = pd.DataFrame(rows)
        t["p_excess_bh"] = stats.false_discovery_control(t["p_excess"].values)
        t["library"] = f"{seg}-segment"
        frames.append(t)
    return pd.concat(frames, ignore_index=True)


def load_table(lib):
    """Mean interaction as a function of the number of further substitutions."""
    frames = []
    for seg, d in lib.items():
        se = SCORE_SE[seg]
        rec = [r for v in contrasts(d["X"], d["y"]).values() for r in v]
        t = pd.DataFrame(rec, columns=["load", "eps"]).groupby("load")["eps"]
        t = t.agg(["mean", "count"]).reset_index()
        t["se"] = 2 * (2 * se) / np.sqrt(t["count"])
        t["library"] = f"{seg}-segment"
        frames.append(t.round(5))
    return pd.concat(frames, ignore_index=True)


# --------------------------------------------------------------------------
# landscape structure
# --------------------------------------------------------------------------
def neighbours(X):
    index = {tuple(row): n for n, row in enumerate(map(tuple, X))}
    nb = []
    for row in X:
        v = []
        for k in range(X.shape[1]):
            h = list(row)
            h[k] = 1 - h[k]
            q = index.get(tuple(h))
            if q is not None:
                v.append(q)
        nb.append(v)
    return nb


def optima_table(lib):
    """Local optima under single-substitution moves and their basins."""
    frames = []
    for seg, d in lib.items():
        X, y = d["X"], d["y"]
        nb = neighbours(X)
        n = len(y)
        dest = np.full(n, -1)
        for start in range(n):
            path, cur = [], start
            while dest[cur] < 0:
                path.append(cur)
                best = max(nb[cur], key=lambda q: y[q]) if nb[cur] else cur
                if y[best] <= y[cur]:
                    dest[cur] = cur
                    break
                cur = best
            end = dest[cur]
            for p in path:
                dest[p] = end
        peaks = [i for i in range(n) if all(y[i] >= y[q] for q in nb[i])]
        counts = pd.Series(dest).value_counts()
        t = pd.DataFrame([dict(peak_index=int(p),
                               score=round(float(y[p]), 5),
                               substitutions=int(X[p].sum()),
                               basin_pct=round(100 * int(counts.get(p, 0)) / n, 3),
                               genotype="".join(map(str, X[p])))
                          for p in peaks]).sort_values("score", ascending=False)
        t["library"] = f"{seg}-segment"
        frames.append(t)
    return pd.concat(frames, ignore_index=True)


def substitution_table(lib):
    """Effect of each substitution on the parental genotype and across backgrounds."""
    frames = []
    for seg, d in lib.items():
        X, y, se = d["X"], d["y"], SCORE_SE[seg]
        index = {tuple(row): n for n, row in enumerate(map(tuple, X))}
        rows = []
        for k, name in enumerate(d["labels"]):
            diffs = []
            for n, row in enumerate(X):
                if row[k] == 0:
                    continue
                h = list(row)
                h[k] = 0
                q = index.get(tuple(h))
                if q is not None:
                    diffs.append(y[n] - y[q])
            diffs = np.array(diffs)
            parent = [0] * X.shape[1]
            on_parent = np.nan
            h = parent.copy()
            h[k] = 1
            if tuple(h) in index and tuple(parent) in index:
                on_parent = float(y[index[tuple(h)]] - y[index[tuple(parent)]])
            rows.append(dict(sub=name,
                             eff_on_parent=round(on_parent, 5),
                             eff_mean_all_bg=round(float(diffs.mean()), 5),
                             eff_min=round(float(diffs.min()), 5),
                             eff_max=round(float(diffs.max()), 5),
                             frac_bg_negative=round(float((diffs < 0).mean()), 5),
                             frac_bg_sig_negative=round(
                                 float((diffs < -2 * np.sqrt(2) * se).mean()), 5),
                             n_bg=len(diffs)))
        t = pd.DataFrame(rows)
        t["library"] = f"{seg}-segment"
        frames.append(t)
    return pd.concat(frames, ignore_index=True)


# --------------------------------------------------------------------------
# models
# --------------------------------------------------------------------------
def design(X, order):
    """Intercept plus one indicator per substitution and per combination."""
    n, L = X.shape
    cols = [np.ones(n)]
    for o in range(1, order + 1):
        for combo in itertools.combinations(range(L), o):
            cols.append(np.prod(X[:, combo], axis=1))
    return np.column_stack(cols)


def cv_r2(Z, y, folds, rng):
    idx = rng.permutation(len(y))
    pred = np.empty(len(y))
    for f in range(folds):
        test = idx[f::folds]
        train = np.setdiff1d(idx, test)
        beta = lstsq(Z[train], y[train], rcond=None)[0]
        pred[test] = Z[test] @ beta
    return 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()


def model_order_table(lib):
    rows = []
    for seg, d in lib.items():
        rng = np.random.default_rng(SEED_CV)
        r2 = {o: cv_r2(design(d["X"], o), d["y"], CV_FOLDS, rng) for o in (1, 2, 3)}
        rows.append(dict(library=f"{seg}-segment",
                         r2_additive=round(r2[1], 4),
                         r2_pairwise=round(r2[2], 4),
                         r2_third_order=round(r2[3], 4),
                         gain_third_over_pairwise=round(r2[3] - r2[2], 4)))
    return pd.DataFrame(rows)


def transfer_table(lib):
    """Pairwise coefficients fitted to low-order data vs to the complete library."""
    rows = []
    for seg, d in lib.items():
        X, y = d["X"], d["y"]
        L = X.shape[1]
        Z = design(X, 2)
        low = X.sum(1) <= 2
        beta_low = lstsq(Z[low], y[low], rcond=None)[0]
        beta_full = lstsq(Z, y, rcond=None)[0]
        start = 1 + L                      # skip intercept and main effects
        a, b = beta_low[start:], beta_full[start:]
        rows.append(dict(library=f"{seg}-segment",
                         n_low_order_genotypes=int(low.sum()),
                         n_pairs=len(a),
                         pearson_r=round(float(np.corrcoef(a, b)[0, 1]), 4)))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# in silico comparison of engineering strategies
# --------------------------------------------------------------------------
def run_strategy(X, y, se, rng, kind, dms=None, n_random=None):
    """One run of one strategy; returns (true score of the pick, measurements)."""
    n, L = X.shape
    observed = y + rng.normal(0, se, n)
    index = {tuple(row): i for i, row in enumerate(map(tuple, X))}
    measured = set()

    if kind == "dms":
        pick = int(np.argmax(X @ dms))
        return float(y[pick]), 1

    if kind == "greedy":
        cur = index[tuple([0] * L)]
        measured.add(cur)
        while True:
            nb = []
            for k in range(L):
                h = list(X[cur])
                h[k] = 1 - h[k]
                q = index.get(tuple(h))
                if q is not None:
                    nb.append(q)
                    measured.add(q)
            best = max(nb, key=lambda q: observed[q])
            if observed[best] > observed[cur]:
                cur = best
            else:
                break
        pool = np.array(sorted(measured))
    else:
        k = X.sum(1)
        if kind == "low2":
            sel = np.where(k <= 2)[0]
        elif kind == "low3":
            sel = np.where(k <= 3)[0]
        elif kind == "random":
            sel = rng.choice(n, n_random, replace=False)
        elif kind == "complete":
            sel = np.arange(n)
        else:
            raise ValueError(kind)
        measured.update(sel.tolist())
        if kind != "complete":
            Z = design(X, 2)
            beta = lstsq(Z[sel], y[sel] + rng.normal(0, se, len(sel)), rcond=None)[0]
            pred = Z @ beta
            mask = np.ones(n, bool)
            mask[sel] = False
            top = np.where(mask)[0][np.argsort(-pred[mask])[:N_VALIDATION]]
            measured.update(top.tolist())
        pool = np.array(sorted(measured))

    # common selection step: re-measure the ten best candidates once
    nominated = pool[np.argsort(-observed[pool])[:N_VALIDATION]]
    second = y[nominated] + rng.normal(0, se, len(nominated))
    pick = nominated[np.argmax((observed[nominated] + second) / 2)]
    return float(y[pick]), len(measured) + len(nominated)


def strategy_table(lib, assignment):
    frames = []
    for seg, d in lib.items():
        X, y, se = d["X"], d["y"], SCORE_SE[seg]
        L = X.shape[1]
        dms = (assignment[assignment.library == seg].sort_values("col")
               ["dms_score"].values.astype(float))
        n_pair = 1 + L + L * (L - 1) // 2
        n_third = n_pair + L * (L - 1) * (L - 2) // 6
        plans = [("combine all beneficial (DMS only)", "dms", None),
                 ("greedy stepwise stacking", "greedy", None),
                 ("low-order design, \u22642 subs.", "low2", None),
                 ("low-order design, \u22643 subs.", "low3", None),
                 ("random design, small n", "random", n_pair),
                 ("random design, larger n", "random", n_third),
                 ("complete library", "complete", None)]
        rng = np.random.default_rng(SEED_STRATEGY)
        rows = []
        for label, kind, size in plans:
            res = [run_strategy(X, y, se, rng, kind, dms=dms, n_random=size)
                   for _ in range(N_STRATEGY_RUNS)]
            vals = np.array([r[0] for r in res])
            meas = np.array([r[1] for r in res])
            rows.append(dict(strategy=label,
                             measurements=int(np.median(meas)),
                             true_score_med=float(np.median(vals)),
                             q10=float(np.percentile(vals, 10)),
                             q90=float(np.percentile(vals, 90)),
                             pct_of_best=100 * float(np.median(vals)) / y.max()))
        t = pd.DataFrame(rows)
        t["library"] = f"{seg}-segment"
        frames.append(t)
    return pd.concat(frames, ignore_index=True)


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", default="data")
    ap.add_argument("--out", default="results")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    lib = load_libraries(args.data)
    assignment = pd.read_csv(os.path.join(args.data, "substitution_assignment.csv"))

    products = {
        "library_summary.csv": library_summary(lib),
        "pairwise_interaction_detectability.csv": pairwise_table(lib),
        "interaction_vs_background_load.csv": load_table(lib),
        "local_optima.csv": optima_table(lib),
        "substitution_background_effects.csv": substitution_table(lib),
        "model_order_cv.csv": model_order_table(lib),
        "coefficient_transfer.csv": transfer_table(lib),
        "strategy_benchmark.csv": strategy_table(lib, assignment),
    }
    for name, table in products.items():
        table.to_csv(os.path.join(args.out, name), index=False)
        print(f"{name}: {len(table)} rows")


if __name__ == "__main__":
    main()
