from __future__ import annotations
from pathlib import Path
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT_DIR = Path(__file__).resolve().parents[1] / "figures" / "local_cohort"
OUT_DIR.mkdir(parents=True, exist_ok=True)
C_R = "#d62728"
C_S = "#2ca02c"
C_I = "#9aa5b1"
C_LOF = "#ff7f0e"
C_REG = "#1f77b4"
C_FS = "#9467bd"
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


def panel_a_isolate_donut(ax) -> None:
    labels = ["Resistant", "Susceptible", "Indeterminate"]
    values = [81, 1043, 341]
    colors = [C_R, C_S, C_I]
    wedges, _ = ax.pie(
        values,
        colors=colors,
        startangle=90,
        counterclock=False,
        wedgeprops=dict(width=0.42, edgecolor="white", linewidth=1.5),
    )
    ax.text(
        0,
        0.06,
        "1,465",
        ha="center",
        va="center",
        fontsize=20,
        fontweight="bold",
        color="#333333",
    )
    ax.text(
        0,
        -0.14,
        "variant carriers",
        ha="center",
        va="center",
        fontsize=9,
        color="dimgrey",
    )
    ax.set_title("A. Isolate-level calls among 1,465 variant carriers")
    legend_items = [
        f"Resistant: n = 81 (5.5%)",
        f"Susceptible: n = 1,043 (71.2%)",
        f"Indeterminate: n = 341 (23.3%)",
    ]
    ax.legend(
        wedges,
        legend_items,
        frameon=False,
        loc="center left",
        bbox_to_anchor=(0.98, 0.5),
        fontsize=9,
    )
    ax.text(
        0,
        -1.25,
        "Resistant isolates account for 1.6% of all 5,000 isolates",
        ha="center",
        va="top",
        fontsize=8.5,
        color="dimgrey",
    )


def panel_b_resistance_evidence(ax) -> None:
    genes = ["$\\it{Rv0678}$", "$\\it{pepQ}$"]
    lof = [48, 22]
    fs = [9, 2]
    y = np.arange(len(genes))[::-1]
    ax.barh(y, lof, height=0.55, color=C_LOF, label="Loss of function")
    ax.barh(y, fs, left=lof, height=0.55, color=C_FS, label="Frameshift")
    for yi, total in zip(y, [57, 24]):
        ax.annotate(
            f"n = {total}",
            (total, yi),
            xytext=(3, 0),
            textcoords="offset points",
            va="center",
            fontsize=9,
        )
    ax.set_yticks(y)
    ax.set_yticklabels(genes, fontsize=10)
    ax.set_xlim(0, 70)
    ax.set_xlabel("Resistance calls (strain–variant pairs, n)")
    ax.set_title("B. Resistance calls by gene and evidence")
    ax.legend(frameon=False, loc="lower right", fontsize=8.5)


def panel_c_structural_lof(ax) -> None:
    genes = ["$\\it{mmpL5}$", "$\\it{mmpS5}$", "$\\it{lpqB}$"]
    values = [24, 3, 2]
    y = np.arange(len(genes))[::-1]
    bars = ax.barh(y, values, height=0.55, color=C_S, alpha=0.85)
    for rect in bars:
        ax.annotate(
            f"n = {int(rect.get_width())}",
            (rect.get_width(), rect.get_y() + rect.get_height() / 2),
            xytext=(3, 0),
            textcoords="offset points",
            va="center",
            fontsize=9,
        )
    ax.set_yticks(y)
    ax.set_yticklabels(genes, fontsize=10)
    ax.set_xlim(0, 30)
    ax.set_xlabel("Loss-of-function variants classified\nas susceptible (n)")
    ax.set_title("C. LoF in structural-integrity genes → susceptible")
    ax.text(
        0.02,
        -0.42,
        "$\\it{mtrA}$/$\\it{mtrB}$: same sensitizing logic. $\\it{glpK}$/$\\it{Rv1979c}$ (exploratory): recorded without a resistance call.",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8,
        color="dimgrey",
    )


def panel_d_variant_outcome(ax) -> None:
    labels = ["Resistant", "Susceptible", "Indeterminate"]
    values = [54, 406, 59]
    colors = [C_R, C_S, C_I]
    y = np.arange(len(labels))[::-1]
    bars = ax.barh(y, values, height=0.55, color=colors)
    for rect in bars:
        ax.annotate(
            f"n = {int(rect.get_width())}",
            (rect.get_width(), rect.get_y() + rect.get_height() / 2),
            xytext=(3, 0),
            textcoords="offset points",
            va="center",
            fontsize=9,
        )
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=10)
    ax.set_xlim(0, 470)
    ax.set_xlabel("Unique variants (n)")
    ax.set_title("D. Variant-level pipeline outcome")
    ax.text(
        0.02,
        -0.42,
        "54 resistance calls = 43 full-pipeline calls + 11 screening-stage frameshifts (9 $\\it{Rv0678}$, 2 $\\it{pepQ}$).",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8,
        color="dimgrey",
    )


def main() -> None:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.6, 4.4))
    panel_a_isolate_donut(ax1)
    panel_b_resistance_evidence(ax2)
    fig.tight_layout()
    stem = "fig_local_cohort_application"
    for ext in ("png", "pdf", "svg"):
        fig.savefig(OUT_DIR / f"{stem}.{ext}")
    plt.close(fig)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.6, 4.2))
    panel_c_structural_lof(ax1)
    panel_d_variant_outcome(ax2)
    fig.tight_layout()
    stem = "figS_local_cohort_supplementary"
    for ext in ("png", "pdf", "svg"):
        fig.savefig(OUT_DIR / f"{stem}.{ext}")
    plt.close(fig)
    print(f"Figures written to: {OUT_DIR}")


if __name__ == "__main__":
    main()
