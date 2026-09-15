import csv
import io
import re
import sys
from collections import Counter, defaultdict
import docx
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image
from docx.oxml.ns import qn

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
RPT = "d:\\python_work\\amr_hunter\\amr-hunter\\reports\\publication"
DOCX = "C:\\Users\\zhang\\Desktop\\深度学习\\3、文章\\正文\\manuscript_EN_02.docx"
OUT_DIR = "d:\\python_work\\amr_hunter\\amr-hunter\\figures"
OUT = OUT_DIR + "\\fig2_combined.png"
fig2 = Image.open(
    OUT_DIR + "\\cryptic_validation\\fig2_cryptic_external_validation.png"
).convert("RGB")
print("fig2 size:", fig2.size, "aspect h/w=%.3f" % (fig2.height / fig2.width))
UP = {"p.Ile67Leu", "p.Phe19Ser", "p.Ser68Asn"}
DN = {"p.Glu55Asp"}
pheno = {}
with open(
    RPT + "\\final_external_bdq_phenotype_validation_isolates.csv", encoding="utf-8-sig"
) as f:
    for r in csv.DictReader(f):
        pheno[r["unique_id"]] = (r["bdq_binary_phenotype"], r.get("bdq_mic") or "")
by = defaultdict(list)
with open(
    RPT + "\\final_external_bdq_phenotype_validation_variant_calls.csv",
    encoding="utf-8-sig",
) as f:
    for r in csv.DictReader(f):
        p = (r.get("amr_hunter_prediction") or "").strip().upper()
        if p not in ("R", "S"):
            p = "UNKNOWN"
        g = (r.get("gene") or "").strip()
        m = (r.get("mutation") or "").strip()
        vc = (r.get("variant_class") or "").strip()
        ev = (r.get("evidence_code") or "").strip()
        if g == "Rv0678" and m in UP:
            p = "R"
        if g == "Rv0678" and m in DN and (p == "R"):
            p = "UNKNOWN"
        if g == "mtrB" and ev == "MTRAB_DEREGULATION":
            p = "S"
        if (
            g == "pepQ"
            and ev == "NEGATIVE_LOF_REVIEWED"
            and (vc == "protein_substitution")
        ):
            p = "UNKNOWN"
        by[r["unique_id"]].append((g, m, p))
iso = {
    u: (
        "R"
        if any((x[2] == "R" for x in rs))
        else "S" if any((x[2] == "S" for x in rs)) else "UNKNOWN"
    )
    for u, rs in by.items()
}
fp = [u for u in pheno if pheno[u][0] == "S" and iso[u] == "R"]


def micnum(s):
    return {"0.5": 0.5, "1.0": 1.0, "2.0": 2.0}.get(
        s, 4.0 if s in (">2", ">1") else 0.0
    )


def miccat(s):
    if s in (">2", ">1"):
        return 4.0
    try:
        return float(s.replace("<=", "").replace(">", ""))
    except ValueError:
        return 0.0


order = [0.008, 0.015, 0.03, 0.06, 0.12, 0.25]
lab = {
    0.008: "≤0.008",
    0.015: "≤0.015",
    0.03: "0.03",
    0.06: "0.06",
    0.12: "0.12",
    0.25: "0.25",
}


def distr(us):
    c = Counter()
    for u in us:
        c[miccat(pheno[u][1])] += 1
    return c


fp_c = distr(fp)
all_s = distr([u for u in pheno if pheno[u][0] == "S"])
fpp = [100 * fp_c.get(k, 0) / len(fp) for k in order]
allp = [100 * all_s.get(k, 0) / sum(all_s.values()) for k in order]


def micbin_res(s):
    s = str(s)
    if s == "0.5":
        return "0.5"
    if s in ("1.0", "1"):
        return "1.0"
    return "≥ 2.0"


r_by = Counter()
t_by = Counter()
for u in pheno:
    if pheno[u][0] == "R":
        b = micbin_res(pheno[u][1])
        r_by[b] += 1
        if iso.get(u) == "R":
            t_by[b] += 1
sens = [
    (b, round(100 * t_by[b] / r_by[b], 1), "%d / %d" % (t_by[b], r_by[b]))
    for b in ("0.5", "1.0", "≥ 2.0")
]
print("G data:", sens)
import numpy as np

TITLE_A = "A. Overall performance against MIC-based phenotypes"
TITLE_F = "F. False-positive MIC shift"


def first_dark_x(arr, y0, y1, x0, x1, thresh=110):
    sub = arr[y0:y1, x0:x1]
    if sub.size == 0:
        return None
    dark = sub.mean(axis=2) < thresh
    rows = np.argwhere(dark.any(axis=1))
    if len(rows) == 0:
        return None
    xs = []
    for r in rows[:60]:
        rr = int(r[0])
        cols = np.argwhere(dark[rr])
        if len(cols):
            xs.append(int(cols[0][0]))
    return min(xs) + x0 if xs else None


plt.rcParams.update(
    {"font.size": 10, "axes.spines.top": False, "axes.spines.right": False}
)
fig = plt.figure(figsize=(12.8, 3.4))
gs = fig.add_gridspec(1, 2, width_ratios=[2.8, 1], wspace=0.8)
ax = fig.add_subplot(gs[0, 0])
y = range(len(order))
b1 = ax.barh(
    [i - 0.19 for i in y],
    fpp,
    height=0.36,
    color="#d6604d",
    alpha=0.92,
    label="False positives (n = 81)",
)
b2 = ax.barh(
    [i + 0.19 for i in y],
    allp,
    height=0.36,
    color="#9fb6cd",
    alpha=0.9,
    label="All susceptible (n = 8,464)",
)
ax.set_yticks(list(y))
ax.set_yticklabels([lab[k] for k in order])
ax.set_ylabel("BDQ MIC (µg/mL)")
ax.set_xlabel("Percentage of isolates (%)")
ax.set_title(TITLE_F)
ax.bar_label(b1, fmt="%.1f", fontsize=8, padding=1)
ax.bar_label(b2, fmt="%.1f", fontsize=8, padding=1)
ax.set_xlim(0, 50)
ax.legend(
    frameon=False,
    fontsize=8.5,
    loc="upper right",
    bbox_to_anchor=(0.96 + 10 / 153.4, 0.9),
)
ax.invert_yaxis()
axF = ax
ax = fig.add_subplot(gs[0, 1])
cuts = [s[0] for s in sens]
vals = [s[1] for s in sens]
ns = [s[2] for s in sens]
bars = ax.bar(cuts, vals, width=0.5, color="#4575b4", alpha=0.92)
for b, v, n in zip(bars, vals, ns):
    ax.text(
        b.get_x() + b.get_width() / 2,
        v + 1.2,
        "%.1f%%\n(%s)" % (v, n),
        ha="center",
        fontsize=9,
        fontweight="bold",
        color="#1a3c6e",
    )
ax.set_ylim(0, 52)
ax.set_ylabel("Sensitivity (%)")
ax.set_xlabel("BDQ MIC (µg/mL)")
ax.set_title("G. Sensitivity by MIC level")
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
axwF = axF.get_window_extent(renderer).width
arrA = np.asarray(fig2.convert("RGB"))
xA = first_dark_x(arrA, 0, 90, 0, arrA.shape[1])
print("A title left px:", xA)
strip_path = OUT_DIR + "\\fig2_mic_fg.png"
anchor = 0.5
gap = int(fig2.width * 0.02)
for it in range(4):
    axF.title.set_position((anchor, 1.0))
    fig.savefig(strip_path, dpi=300, bbox_inches="tight")
    strip = Image.open(strip_path).convert("RGB")
    W = fig2.width
    k = W / strip.width
    strip_r = strip.resize((W, int(strip.height * k)), Image.LANCZOS)
    canvas = Image.new("RGB", (W, fig2.height + gap + strip_r.height), "white")
    canvas.paste(fig2, (0, 0))
    canvas.paste(strip_r, (0, fig2.height + gap))
    arrC = np.asarray(canvas)
    y0 = fig2.height + gap + 6
    xF = first_dark_x(arrC, y0, y0 + 90, 0, 1200)
    print("iter", it, "xF:", xF, "xA:", xA, "anchor: %.4f" % anchor)
    if xF is None or xA is None or abs(xF - xA) <= 2:
        canvas.save(OUT)
        break
    anchor += (xA - xF) / (axwF * 3 * k)
plt.close(fig)
print(
    "combined saved:",
    OUT,
    canvas.size,
    "aspect h/w=%.3f" % (canvas.height / canvas.width),
)
