"""New Fig 3 - 4 columns x 3 rows; the right half of row 2 is split into two
stacked horizontal panels for the strategy comparison.
  row 1  a,b  editing score vs number of substitutions (N, C)
         c,d  local optima (N, C)
  row 2  e,f  mean pairwise interaction vs background load (N, C)
         g    strategies, N (horizontal, upper)
         h    strategies, C (horizontal, lower)
  row 3  i,j  per-substitution effects (N, C)
Requires in-kernel: LIB, effs, allbasin, BEN, LOAD, SUBS,
apply_figure_style, panel_letter, META_GREY.
"""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt

CB = {"N": "#1f5fa9", "C": "#d1913c", "seq": "#8c8c8c", "dms": "#b3452c",
      "pick": "#000000",
      "low": "#d1913c", "full": "#1f5fa9", "drop": "#b0b0b0"}
NICE = {"N": "N-segment", "C": "C-segment"}
SHOW_LETTERS = False   # panel letters are added by hand at layout
RHO = {"N": (0.45, "0.17"), "C": (0.82, "0.004")}
FS = 6.8  # every in-panel text: annotations, legends, category labels
PARENT = {"N": "wild type", "C": "Cas12j-N7"}
CARRIED = {"N": "Cas12j-N7", "C": "evoCas12j"}
STRAT = [("combine all beneficial (DMS only)", "DMS only", "D", CB["dms"]),
         ("greedy stepwise stacking", "greedy", "^", CB["seq"]),
         ("low-order design, \u22642 subs.", "low-order \u22642", "s", CB["low"]),
         ("low-order design, \u22643 subs.", "low-order \u22643", "s", CB["low"]),
         ("random design, small n", "random, small", "o", CB["full"]),
         ("random design, larger n", "random, large", "o", CB["full"]),
         ("complete library", "complete", "*", CB["full"])]


def build():
    plt.close("all")
    apply_figure_style(sizes=(8.0, 6.8, 6.8))
    fig = plt.figure(figsize=(7.2, 7.2))
    gs = fig.add_gridspec(3, 4, hspace=0.60, wspace=0.55, height_ratios=[1.0, 1.0, 1.05])
    left_axes = []
    seq = iter("abcdefghij")
    def letter(ax):
        if SHOW_LETTERS:
            panel_letter(ax, next(seq), case="upper")

    # row 1, cols 0-1: score vs number of substitutions
    ax_ab = None
    for i, nm in enumerate(("N", "C")):
        X, y, se = LIB[nm]
        ax = fig.add_subplot(gs[0, i], sharey=ax_ab); ax_ab = ax_ab or ax; k = X.sum(1)
        if i == 1:
            _p = ax.get_position(); ax.set_position([_p.x0 - 0.045, _p.y0, _p.width, _p.height])
        ks = [v for v in range(X.shape[1] + 1) if (k == v).sum() > 0]
        data = [y[k == v] for v in ks]
        for b in ax.violinplot(data, positions=ks, widths=0.9, showextrema=False)['bodies']:
            b.set_facecolor(CB[nm]); b.set_edgecolor("none"); b.set_alpha(0.32)
        ax.plot(ks, [d.max() for d in data], "o-", color=CB[nm], ms=2.5, lw=1)
        bi = int(np.argmax(y)); ax.plot([k[bi]], [y[bi]], "*", color=CB["pick"], ms=9, zorder=5)
        ax.annotate(CARRIED[nm], (k[bi], y[bi]), xytext=(0.1, y[bi] * (1.09 if nm == "N" else 1.75)),
                    color=CB["pick"], fontsize=FS, ha="left", va="center",
                    arrowprops=dict(arrowstyle="-", color=CB["pick"], lw=0.6,
                                    shrinkA=2, shrinkB=5))
        if i == 0:
            left_axes.append(ax)
            ax.set_xlabel("Substitutions carried"); ax.xaxis.set_label_coords(1.07, -0.18)
        else:
            ax.set_xlabel("")
        ax.set_ylim(top=max(LIB[k2][1].max() for k2 in LIB) * 1.17)
        if i == 0:
            ax.set_ylabel("Editing score", labelpad=1)
        else:
            ax.tick_params(labelleft=False)
        if i == 0:
            letter(ax)

    # row 1, cols 2-3: local optima
    ax_cd = None
    for i, nm in enumerate(("N", "C")):
        X, y, se = LIB[nm]
        ax = fig.add_subplot(gs[0, 2 + i], sharey=ax_cd); ax_cd = ax_cd or ax
        _p = ax.get_position()
        ax.set_position([_p.x0 + (0.012 if i == 0 else -0.033), _p.y0, _p.width, _p.height])
        B = allbasin[nm]; bi = int(np.argmax(y))
        ax.plot(B["score"], B["basin_pct"], "o", mfc="none", mec=CB[nm], ms=3, mew=0.8, zorder=3)
        row = B[B["peak"] == bi]
        ax.plot(row["score"], row["basin_pct"], "D", color=CB["pick"], ms=4.5, zorder=5)
        if nm == "C":
            w = B.iloc[0]
            ax.plot([w["score"]], [w["basin_pct"]], "o", mfc="white", mec="k", ms=5.5, mew=1.0, zorder=5)
            ax.annotate("widest basin", (w["score"], w["basin_pct"]),
                        xytext=(w["score"] - 0.030, w["basin_pct"] * 1.6), fontsize=FS)
        ax.set_yscale("log"); ax.set_ylim(0.4, 95)
        if i == 0:
            ax.set_xlabel("Peak editing score"); ax.xaxis.set_label_coords(1.07, -0.18)
        else:
            ax.set_xlabel("")
        if i == 0:
            ax.set_ylabel("Basin (% of library)", labelpad=1)
        else:
            ax.tick_params(labelleft=False)
        if i == 0:
            letter(ax)

    # row 2, cols 0-1: parental-background cycle vs mean over backgrounds
    ax_ef = None
    for i, nm in enumerate(("N", "C")):
        X, y, se = LIB[nm]
        ax = fig.add_subplot(gs[1, i], sharey=ax_ef); ax_ef = ax_ef or ax
        if i == 1:
            _p = ax.get_position(); ax.set_position([_p.x0 - 0.045, _p.y0, _p.width, _p.height])
        g = PAIR[PAIR["library"] == f"{nm}-segment"].copy()
        g["se_mean"] = 2 * se / np.sqrt(g["n_bg"])
        thr = 2 * (2 * se)
        lim = max(g["eps_wt"].abs().max(), g["eps_mean"].abs().max()) * 1.3
        for xv in (-thr, thr):
            ax.axvline(xv, color=META_GREY, lw=0.7, ls=":", zorder=1)
        ax.plot([-lim, lim], [-lim, lim], color="k", lw=0.6, ls="--", zorder=1)
        ax.axhline(0, color="k", lw=0.6); ax.axvline(0, color="k", lw=0.6)
        sig = g["sig_complete"].values
        ax.errorbar(g["eps_wt"][~sig], g["eps_mean"][~sig], yerr=2 * g["se_mean"][~sig],
                    fmt="o", ms=2.6, mfc="none", mec=CB["drop"], ecolor=CB["drop"],
                    elinewidth=0.5, lw=0, zorder=2)
        ax.errorbar(g["eps_wt"][sig], g["eps_mean"][sig], yerr=2 * g["se_mean"][sig],
                    fmt="o", ms=3.0, color=CB[nm], ecolor=CB[nm], elinewidth=0.5, lw=0, zorder=3)
        ax.set_xlim(-lim, lim)
        if i == 0:
            ax.set_ylim(-0.22, 0.22)
            ax.set_ylabel("\u03b5, mean over all backgrounds", labelpad=1)
        else:
            ax.tick_params(labelleft=False)
        if i == 0:
            left_axes.append(ax)
            ax.set_xlabel("\u03b5, parental background"); ax.xaxis.set_label_coords(1.07, -0.20)
        else:
            ax.set_xlabel("")
        _first_ef = (i == 0)
        ax.text(0.03, 0.96, f"{int(sig.sum())}/{len(g)} resolved",
                transform=ax.transAxes, fontsize=FS, color=CB[nm], va="top")
        if i == 0:
            letter(ax)

    # row 2, cols 2-3: strategies, both libraries merged into one panel
    ax = fig.add_subplot(gs[1, 2:4])
    _p = ax.get_position(); ax.set_position([_p.x0 + 0.075, _p.y0, _p.width - 0.075, _p.height])
    ys = np.arange(len(STRAT))
    ax.axvline(100, color=META_GREY, lw=0.8, ls=":")
    labs = []
    for yy, (key, short, mk, col) in enumerate(zip(*zip(*[(k, s_, m_, c_) for k, s_, m_, c_ in STRAT]))) if False else zip(ys, STRAT):
        key, short, mk, _ = yy if False else (STRAT[int(yy)][0], STRAT[int(yy)][1], STRAT[int(yy)][2], None)
        ns = []
        for nm, dy in (("N", -0.19), ("C", 0.19)):
            X, y, se = LIB[nm]
            r = BEN[BEN["library"] == f"{nm}-segment"].set_index("strategy").loc[key]
            best = y.max()
            m, q1, q9 = (100 * r[c] / best for c in ("true_score_med", "q10", "q90"))
            ax.plot([q1, q9], [yy + dy] * 2, color=CB[nm], lw=1.2, solid_capstyle="round", zorder=2)
            ax.plot([m], [yy + dy], mk, color=CB[nm], ms=9 if mk == "*" else 4.6, zorder=3)
            ns.append(int(r["measurements"]))
        labs.append(short); ax.text(113, yy, f"{ns[0]:,} / {ns[1]:,}", fontsize=FS,
                                    color="#555555", va="center", ha="left")
    ax.set_yticks(ys); ax.set_yticklabels(labs, fontsize=FS)
    ax.set_xticks([40, 60, 80, 100])
    ax.invert_yaxis(); ax.set_xlim(34, 142); ax.set_ylim(len(STRAT) - 0.4, -0.6)
    ax.set_xlabel("Activity of the selected variant\n(% of best measured genotype)")
    ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
    letter(ax)

    # row 3: per-substitution effects
    for i, nm in enumerate(("N", "C")):
        X, y, se = LIB[nm]
        ax = fig.add_subplot(gs[2, 2 * i:2 * i + 2])
        if i == 1:
            _q2 = ax.get_position()
            ax.set_position([_q2.x0 + 0.022, _q2.y0, _q2.width - 0.022, _q2.height])
        _q = ax.get_position()
        if nm == "N":
            ax.set_position([_q.x0, _q.y0, _q.width - 0.035, _q.height])
        else:
            ax.set_position([_q.x0 + 0.040, _q.y0, _q.width - 0.040, _q.height])
        S = SUBS[SUBS["library"] == f"{nm}-segment"].copy().reset_index(drop=True)
        S["j"] = range(len(S)); S = S.sort_values("eff_mean_all_bg")
        adopted = set(np.where(X[int(np.argmax(y))] == 1)[0]); pos = np.arange(len(S))
        vp = ax.violinplot([effs[nm][j] for j in S["j"]], positions=pos, widths=0.8,
                           showextrema=False)
        for bb, j in zip(vp['bodies'], S["j"]):
            bb.set_facecolor(CB[nm] if j in adopted else CB["drop"])
            bb.set_alpha(0.5); bb.set_edgecolor("none")
        ax.plot(pos, S["eff_mean_all_bg"], "o", color="k", ms=3, zorder=4,
                label="mean over all backgrounds")
        ax.plot(pos, S["eff_on_parent"], "D", mfc="white", mec="#000000", ms=4.5, zorder=4,
                label=f"effect on the parent ({PARENT[nm]})")
        ax.axhline(0, color="k", lw=0.7)
        ax.set_xticks(pos); ax.set_xticklabels(S["sub"], rotation=45, ha="right", fontsize=FS)
        for t, j in zip(ax.get_xticklabels(), S["j"]):
            t.set_color(CB[nm] if j in adopted else "#6b6b6b")
        if i == 0:
            left_axes.append(ax)
        rr, pp = RHO[nm]
        ax.text(0.98, 0.04, f"Spearman \u03c1 = {rr:.2f}, P = {pp}", transform=ax.transAxes,
                fontsize=FS, ha="right", va="bottom", color="#444444")
        ax.set_ylabel("Change in editing score", labelpad=1)
        ax.set_ylim(-0.45, 0.60) if nm == "N" else ax.set_ylim(-0.17, 0.25)
        ax.legend(frameon=False, fontsize=FS, loc="upper left", handletextpad=0.4, borderpad=0.1)
        letter(ax)

    W = 0.016
    for il, ir in ((0, 1), (4, 5)):
        pl, pr = fig.axes[il].get_position(), fig.axes[ir].get_position()
        fig.axes[il].set_position([pl.x0, pl.y0, pl.width + W, pl.height])
        fig.axes[ir].set_position([pr.x0 + W, pr.y0, pr.width + W, pr.height])

    x0 = min(a.get_position().x0 for a in left_axes)
    for a in left_axes:
        p = a.get_position(); a.set_position([x0, p.y0, p.x1 - x0, p.height])

    fig.savefig("fig3_new.png", dpi=400, bbox_inches="tight")
    fig.savefig("fig3_new.pdf", bbox_inches="tight")
    return fig


fig3 = build()
