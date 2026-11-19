"""Figure S3 - pairwise interactions in both combinatorial libraries.
  a, b  N-segment: interaction matrix on the parental background / averaged over backgrounds
  c, d  C-segment: the same two matrices
  e     mean interaction vs number of other substitutions in the background
  f     pairwise coefficients fitted to <=2-substitution data vs the complete library
Requires in-kernel: LIB, PAIR, LOAD, coef, labN, labC, apply_figure_style, panel_letter, META_GREY.
"""
import numpy as np, matplotlib as mpl, matplotlib.pyplot as plt

CB = {"N": "#1f5fa9", "C": "#d1913c"}
NICE = {"N": "N-segment", "C": "C-segment"}
LABS = {"N": labN, "C": labC}
SHOW_LETTERS = False   # panel letters are added by hand at layout
FS = 6.8


def _mat(nm, col):
    labs = LABS[nm]; L = len(labs)
    M = np.full((L, L), np.nan)
    for _, r in PAIR[PAIR["library"] == f"{nm}-segment"].iterrows():
        i, j = (int(v) - 1 for v in str(r["pair"]).split("-"))
        M[i, j] = M[j, i] = r[col]
    return M


def build():
    plt.close("all")
    apply_figure_style(sizes=(8.0, 6.8, 6.8))
    fig = plt.figure(figsize=(7.2, 8.6))
    gs = fig.add_gridspec(3, 2, hspace=0.45, wspace=0.45, height_ratios=[1.25, 1.25, 1.0])
    seq = iter("abcdef")
    def letter(ax):
        if SHOW_LETTERS:
            panel_letter(ax, next(seq), case="upper")

    for row, nm in enumerate(("N", "C")):
        labs = LABS[nm]
        v = max(np.nanmax(np.abs(_mat(nm, "eps_wt"))), np.nanmax(np.abs(_mat(nm, "eps_mean"))))
        for col, (key, what) in enumerate((("eps_wt", "parental background"),
                                           ("eps_mean", "mean over backgrounds"))):
            ax = fig.add_subplot(gs[row, col])
            M = _mat(nm, key)
            im = ax.imshow(M, cmap="RdBu_r", vmin=-v, vmax=v)
            ax.set_xticks(range(len(labs))); ax.set_yticks(range(len(labs)))
            ax.set_xticklabels(labs, rotation=90, fontsize=FS)
            ax.set_yticklabels(labs, fontsize=FS)
            ax.tick_params(length=0)
            for s in ax.spines.values():
                s.set_visible(False)
            cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
            cb.ax.tick_params(labelsize=FS, length=2)
            cb.set_label("Interaction \u03b5", fontsize=FS, labelpad=2)
            letter(ax)

    # e: interaction vs background load
    ax = fig.add_subplot(gs[2, 0])
    for nm in ("N", "C"):
        D = LOAD[LOAD["library"] == f"{nm}-segment"].sort_values("load")
        ax.errorbar(D["load"], D["mean"], yerr=2 * D["se"], fmt="o-", color=CB[nm],
                    ms=3, lw=1, capsize=1.5, label=NICE[nm])
    ax.axhline(0, color="k", lw=0.7)
    ax.set_xlabel("Other substitutions in the background")
    ax.set_ylabel("Mean interaction \u03b5", labelpad=1)
    ax.legend(frameon=False, fontsize=FS, loc="lower right", handletextpad=0.4, borderpad=0.1)
    letter(ax)

    # f: coefficient transfer
    ax = fig.add_subplot(gs[2, 1])
    _p = ax.get_position(); ax.set_position([_p.x0 + 0.02, _p.y0, _p.width - 0.02, _p.height])
    lim = 0
    for nm in ("N", "C"):
        lo, fu, r, nlow = coef[nm]
        lim = max(lim, np.abs(lo).max(), np.abs(fu).max())
        ax.plot(lo, fu, "o", mfc="none", mec=CB[nm], ms=3, mew=0.8,
                label=f"{NICE[nm]}: r = {r:.2f}")
    lim *= 1.1
    ax.plot([-lim, lim], [-lim, lim], color="k", lw=0.6, ls="--", zorder=1)
    ax.axhline(0, color="k", lw=0.6); ax.axvline(0, color="k", lw=0.6)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("Coefficient, \u22642-substitution data")
    ax.set_ylabel("Coefficient, complete library", labelpad=1)
    ax.legend(frameon=False, fontsize=FS, loc="upper left", handletextpad=0.4, borderpad=0.1)
    letter(ax)

    fig.savefig("figS3_new.png", dpi=400, bbox_inches="tight")
    fig.savefig("figS3_new.pdf", bbox_inches="tight")
    return fig


figS3 = build()
