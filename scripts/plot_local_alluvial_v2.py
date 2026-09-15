from __future__ import annotations
import csv
import sys
from collections import Counter
from pathlib import Path
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, PathPatch
from matplotlib.path import Path as MplPath

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
BASE = Path("C:\\Users\\zhang\\Desktop\\深度学习\\3、文章\\内部数据")
ANA = BASE / "outdir" / "analysis"
OUT_DIR = Path(__file__).resolve().parents[1] / "figures" / "local_cohort"
C_R, C_S, C_I, C_P = ("#d62728", "#2ca02c", "#9aa5b1", "#dfe3e8")
plt.rcParams.update(
    {"font.size": 10, "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight"}
)
vcalls = Counter()
with (ANA / "variant_calls_new.csv").open(encoding="utf-8") as f:
    for row in csv.DictReader(f):
        if row["gene"] not in {
            "atpE",
            "Rv0678",
            "mmpS5",
            "mmpL5",
            "pepQ",
            "lpqB",
            "mtrA",
            "mtrB",
            "Rv1979c",
            "glpK",
        }:
            continue
        t = row["variant_type"]
        if not (
            t.startswith("Nonsynonymous")
            or "frameshift" in t
            or t.startswith("Intergenic")
        ):
            continue
        vcalls[row["call"]] += 1
print("variant calls:", dict(vcalls))
n_var = sum(vcalls.values())
r, s, i, p = (vcalls["R"], vcalls["S"], vcalls["INDET"], vcalls["PENDING"])
icalls = Counter()
with (ANA / "isolate_calls_new.csv").open(encoding="utf-8") as f:
    for row in csv.DictReader(f):
        if row["if_phylogeny"] == "1":
            icalls[row["call"]] += 1
print("isolate calls (pure):", dict(icalls))
n_iso = sum(icalls.values())
ir, is_, ii = (icalls["R"], icalls["S"], icalls["INDET"])


def flow(ax, x0, w0, yt0, yb0, x1, w1, yt1, yb1, color, alpha=0.55):
    xm = (x0 + w0 + x1) / 2
    verts = [
        (x0 + w0, yt0),
        (xm, yt0),
        (x1, yt1),
        (x1, yb1),
        (xm, yb0),
        (x0 + w0, yb0),
        (0, 0),
    ]
    codes = [
        MplPath.MOVETO,
        MplPath.CURVE3,
        MplPath.CURVE3,
        MplPath.LINETO,
        MplPath.CURVE3,
        MplPath.CURVE3,
        MplPath.CLOSEPOLY,
    ]
    ax.add_patch(
        PathPatch(
            MplPath(verts, codes),
            facecolor=color,
            alpha=alpha,
            edgecolor="none",
            zorder=2,
        )
    )


def segment(ax, x, w, ytop, ybot, color, name, value, pct, hatch=None):
    ax.add_patch(
        FancyBboxPatch(
            (x, ybot),
            w,
            ytop - ybot,
            boxstyle="round,pad=0,rounding_size=0.004",
            facecolor=color,
            edgecolor="white",
            linewidth=1.2,
            hatch=hatch,
            zorder=3,
        )
    )
    ymid = (ytop + ybot) / 2
    ax.text(
        x + w / 2,
        ytop - 0.04,
        name,
        ha="center",
        va="center",
        fontsize=8.5,
        color="white",
        fontweight="bold",
        zorder=5,
    )
    vlabel = value if isinstance(value, str) else f"{value:,}"
    plabel = "—" if pct is None else f"({pct:.1f}%)"
    ax.text(
        x + w / 2,
        ymid + 0.015,
        vlabel,
        ha="center",
        va="center",
        fontsize=11,
        fontweight="bold",
        color="white",
        zorder=5,
    )
    ax.text(
        x + w / 2,
        ymid - 0.042,
        plabel,
        ha="center",
        va="center",
        fontsize=7.5,
        color="white",
        alpha=0.95,
        zorder=5,
    )


fig, ax = plt.subplots(figsize=(10.6, 7.0))
ax.set_xlim(-0.06, 1.02)
ax.set_ylim(0.0, 1.16)
ax.axis("off")
x0, x1, wcol = (0.04, 0.65, 0.14)
y = 1.0
seg1 = []
for v, c, name in (
    (r, C_R, "Resistant"),
    (s, C_S, "Susceptible"),
    (i, C_I, "Indeterminate"),
    (p, C_P, "Pending"),
):
    ytop, ybot = (y, y - v / n_var)
    seg1.append((ytop, ybot, c, name, v))
    y = ybot
y = 1.0
seg2 = []
for v, c, name in (
    (ir, C_R, "Resistant"),
    (is_, C_S, "Susceptible"),
    (ii, C_I, "Indeterminate"),
):
    ytop, ybot = (y, y - v / n_iso)
    seg2.append((ytop, ybot, c, name, v))
    y = ybot
for (vt, vb, c, nm, v), (it, ib, c2, nm2, iv) in zip(seg1[:3], seg2[:3]):
    flow(ax, x0, wcol, vt, vb, x1, wcol, it, ib, c, 0.55)
PH = {"Susceptible", "Indeterminate"}
for vt, vb, c, nm, v in seg1:
    segment(
        ax,
        x0,
        wcol,
        vt,
        vb,
        c,
        nm,
        "XXX" if nm in PH else v,
        None if nm in PH else 100 * v / n_var,
    )
for it, ib, c, nm, v in seg2:
    segment(
        ax,
        x1,
        wcol,
        it,
        ib,
        c,
        nm,
        "XXX" if nm in PH else v,
        None if nm in PH else 100 * v / n_iso,
    )
ax.text(
    x0 + wcol / 2,
    1.055,
    "Variant-level calls",
    ha="center",
    va="bottom",
    fontsize=11,
    fontweight="bold",
)
ax.text(
    x0 + wcol / 2,
    1.012,
    f"(n = {n_var:,})",
    ha="center",
    va="bottom",
    fontsize=8,
    color="dimgrey",
)
ax.text(
    x1 + wcol / 2,
    1.055,
    "Isolate-level calls",
    ha="center",
    va="bottom",
    fontsize=11,
    fontweight="bold",
)
ax.text(
    x1 + wcol / 2,
    1.012,
    f"(n = {n_iso:,})",
    ha="center",
    va="bottom",
    fontsize=8,
    color="dimgrey",
)
out = OUT_DIR / "fig_alluvial_application"
for ext in ("png", "pdf", "svg"):
    fig.savefig(f"{out}.{ext}")
plt.close(fig)
print(f"saved {out}.png/.pdf/.svg")
print(f"flows: R={r} S={s} INDET={i} PENDING={p} | isolates R={ir} S={is_} INDET={ii}")
