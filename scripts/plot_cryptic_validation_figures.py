from __future__ import annotations
from pathlib import Path
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT_DIR = Path(__file__).resolve().parents[1] / "figures" / "cryptic_validation"
OUT_DIR.mkdir(parents=True, exist_ok=True)
C_RECALL = "#1f77b4"
C_PRECISION = "#d62728"
C_CV = "#2ca02c"
C_HOLDOUT = "#ff7f0e"
C_BAR = "#4c72b0"
A_SHIFT_PX = 236.2
B_SHIFT_PX = 118.1
plt.rcParams.update(
    {
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
    }
)


def panel_a_overall_performance(ax) -> None:
    metrics = ["Recall", "Specificity", "Precision"]
    values = [33.8, 99.04, 22.9]
    colors = [C_RECALL, "#2ca02c", C_PRECISION]
    bars = ax.bar(metrics, values, color=colors, width=0.55)
    labels = ["33.8%", "99.04%", "22.9%"]
    for rect, lab in zip(bars, labels):
        ax.annotate(
            lab,
            (rect.get_x() + rect.get_width() / 2, rect.get_height()),
            xytext=(0, 2),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=9,
        )
    ax.annotate(
        "95% CI 23.1–46.2",
        (0, 33.8),
        xytext=(0, 20),
        textcoords="offset points",
        ha="center",
        fontsize=8,
        color="dimgrey",
    )
    ax.set_ylim(0, 108)
    ax.set_ylabel("Percent")
    ax.set_title("A. Overall performance against MIC-based phenotypes")


def panel_b_fn_mic(ax) -> None:
    bins = ["0.5", "1.0", "≥2.0"]
    counts = [28, 9, 10]
    pct = [59.6, 19.1, 21.3]
    bars = ax.bar(bins, counts, color=C_BAR, width=0.62, alpha=0.9)
    for rect, p in zip(bars, pct):
        ax.annotate(
            f"n = {int(rect.get_height())}\n({p:.1f}%)",
            (rect.get_x() + rect.get_width() / 2, rect.get_height()),
            xytext=(0, 2),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=9,
        )
    ax.set_ylim(0, 35)
    ax.set_ylabel("False-negative isolates (n)")
    ax.set_xlabel("BDQ MIC (µg/mL)")
    ax.set_title("B. MICs of the 47 false-negative isolates")


def panel_c_fp_composition(ax) -> None:
    labels = [
        "$\\it{Rv0678}$:p.Asp47fs",
        "$\\it{Rv0678}$:p.Met146Thr",
        "$\\it{Rv0678}$:p.Cys46fs",
        "$\\it{Rv0678}$:p.Arg90Cys",
        "$\\it{Rv0678}$:p.Leu117Arg",
        "$\\it{mtrB}$ & $\\it{pepQ}$",
        "Other",
    ]
    counts = [15, 11, 7, 6, 3, 4, 35]
    colors = [C_BAR] * 6 + ["#9aa5b1"]
    y = np.arange(len(labels))[::-1]
    bars = ax.barh(y, counts, color=colors, height=0.6)
    for rect in bars:
        ax.annotate(
            f"n = {int(rect.get_width())}",
            (rect.get_width(), rect.get_y() + rect.get_height() / 2),
            xytext=(3, 0),
            textcoords="offset points",
            va="center",
            fontsize=8.5,
        )
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.set_xlim(0, 40)
    ax.set_xlabel("False-positive isolates (n)")
    ax.set_title("C. Composition of the 81 false positives")


def combined_figure() -> None:
    fig = plt.figure(figsize=(12.8, 7.6))
    gs = fig.add_gridspec(2, 3, hspace=0.65, wspace=0.8)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2])
    ax_d = fig.add_subplot(gs[1, 0:2])
    ax_e = fig.add_subplot(gs[1, 2])
    panel_a_overall_performance(ax_a)
    panel_b_fn_mic(ax_b)
    panel_c_fp_composition(ax_c)
    panel_a_mtrb_enrichment(ax_d)
    panel_b_supervised_ceiling(ax_e)
    if A_SHIFT_PX or B_SHIFT_PX:
        fig_width_px = fig.get_size_inches()[0] * fig.dpi
        for shift_px, ax in ((A_SHIFT_PX, ax_a), (B_SHIFT_PX, ax_b)):
            if shift_px:
                pos = ax.get_position()
                ax.set_position(
                    [pos.x0 - shift_px / fig_width_px, pos.y0, pos.width, pos.height]
                )
    ax_d.set_title("D. $\\it{mtrB}$ combination enrichment")
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    a_bb = ax_a.title.get_window_extent(renderer)
    d_bb = ax_d.title.get_window_extent(renderer)
    axes_width_px = ax_d.get_window_extent(renderer).width
    dx_frac = (a_bb.x0 - d_bb.x0) / axes_width_px
    ax_d.title.set_position((0.5 + dx_frac, 1.0))
    stem = "fig2_cryptic_external_validation"
    for ext in ("png", "pdf", "svg"):
        fig.savefig(OUT_DIR / f"{stem}.{ext}", bbox_inches="tight")
    plt.close(fig)


def panel_a_mtrb_enrichment(ax) -> None:
    rows = [
        ("$\\it{mtrB}$:p.Met517Leu alone", 0.78, None, "P = 0.56"),
        (
            "$\\it{mtrB}$:p.Met517Leu\n+ secondary variants",
            3.95,
            (2.43, 6.44),
            "P = 3.3 × 10⁻⁸",
        ),
        ("$\\it{mtrB}$:p.Met517Leu + $\\it{glpK}$", 2.46, None, "P = 1.0 × 10⁻³"),
        ("$\\it{mtrB}$:p.Met517Leu + $\\it{Rv1979c}$", 2.5, None, "P = 1.1 × 10⁻³"),
        ("$\\it{mtrB}$:p.Pro18Ser\n(negative control)", 0.58, None, "P = 0.06"),
    ]
    y = np.arange(len(rows))[::-1]
    ax.axvline(1.0, color="grey", ls="--", lw=1, zorder=0)
    for yi, (label, orv, ci, ptext) in zip(y, rows):
        ax.scatter(orv, yi, color=C_BAR, s=34, zorder=3)
        if ci is not None:
            ax.plot(ci, [yi, yi], color=C_BAR, lw=1.6, zorder=2)
            ax.annotate(
                f"{orv:.2f} ({ci[0]:.2f}–{ci[1]:.2f})",
                (orv, yi),
                xytext=(0, 12),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=8.5,
            )
        else:
            ax.annotate(
                f"{orv:.2f}",
                (orv, yi),
                xytext=(6, 0),
                textcoords="offset points",
                va="center",
                fontsize=8.5,
            )
        if ptext:
            ax.annotate(
                ptext, (9.3, yi), ha="right", va="center", fontsize=8, color="dimgrey"
            )
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=8.5)
    ax.set_xscale("log")
    ax.set_xlim(0.2, 9.5)
    ax.set_xticks([0.3, 1.0, 3.0, 8.0])
    ax.set_xticklabels(["0.3", "1", "3", "8"])
    ax.set_xlabel("Odds ratio for phenotypic BDQ resistance (log scale)")


def panel_b_supervised_ceiling(ax) -> None:
    labels = ["Specificity\n≥ 99%", "Specificity\n≥ 95%", "Specificity\n≥ 90%"]
    cv = [26.5, 41.1, 46.5]
    holdout = [23.9, 36.6, 49.3]
    x = np.arange(len(labels))
    w = 0.34
    b1 = ax.bar(x - w / 2, cv, w, color=C_CV, label="Cross-validation")
    b2 = ax.bar(x + w / 2, holdout, w, color=C_HOLDOUT, label="Site-holdout")
    for bars in (b1, b2):
        for rect in bars:
            ax.annotate(
                f"{rect.get_height():.1f}%",
                (rect.get_x() + rect.get_width() / 2, rect.get_height()),
                xytext=(0, 2),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=9,
            )
    ax.axhline(49.3, color="dimgrey", ls="--", lw=1)
    ax.annotate(
        "≈ 49.3% ceiling",
        (2.48, 51.0),
        xycoords="data",
        ha="right",
        va="bottom",
        fontsize=8.5,
        color="dimgrey",
    )
    ax.axhline(33.8, color=C_RECALL, ls=":", lw=1.2)
    ax.annotate(
        "AMR-Hunter zero-shot recall (33.8%)",
        (0.06, 35.4),
        xycoords="data",
        ha="left",
        va="bottom",
        fontsize=8.5,
        color=C_RECALL,
    )
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8.5)
    ax.set_ylim(0, 64)
    ax.set_ylabel("BDQ-resistant recall (%)")
    ax.set_title("E. Supervised genotype-only ceiling")
    ax.legend(frameon=False, loc="upper left", fontsize=9)


def main() -> None:
    combined_figure()
    print(f"Figures written to: {OUT_DIR}")
    for p in sorted(OUT_DIR.iterdir()):
        if p.is_file():
            print(f"  {p.name} ({p.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
