from __future__ import annotations
import csv
import math
from collections import defaultdict
from pathlib import Path
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import binomtest, chi2

RPT = Path("d:\\python_work\\amr_hunter\\amr-hunter\\reports\\publication")
CALLS = "C:\\Users\\zhang\\Desktop\\深度学习\\3、文章\\内部数据\\bdq_calls.csv"
OUT_DIR = Path(__file__).resolve().parents[1] / "figures" / "method_comparison"
OUT_DIR.mkdir(parents=True, exist_ok=True)
UPGRADE = {"p.Ile67Leu", "p.Phe19Ser", "p.Ser68Asn"}
DOWNGRADE = {"p.Glu55Asp"}


def pred(v: str) -> str:
    v = (v or "").strip().upper()
    if v in {"R", "RESISTANT"}:
        return "R"
    if v in {"S", "SENSITIVE", "SUSCEPTIBLE"}:
        return "S"
    return "UNKNOWN"


by_iso: dict[str, list[str]] = defaultdict(list)
with (RPT / "final_external_bdq_phenotype_validation_variant_calls.csv").open(
    newline="", encoding="utf-8-sig"
) as f:
    for r in csv.DictReader(f):
        p = pred(r.get("amr_hunter_prediction"))
        g = (r.get("gene") or "").strip()
        m = (r.get("mutation") or "").strip()
        vc = (r.get("variant_class") or "").strip()
        ev = (r.get("evidence_code") or "").strip()
        if g == "Rv0678" and m in UPGRADE:
            p = "R"
        if g == "Rv0678" and m in DOWNGRADE and (p == "R"):
            p = "UNKNOWN"
        if g == "mtrB" and ev == "MTRAB_DEREGULATION":
            p = "S"
        if (
            g == "pepQ"
            and ev == "NEGATIVE_LOF_REVIEWED"
            and (vc == "protein_substitution")
        ):
            p = "UNKNOWN"
        by_iso[r["unique_id"]].append(p)
amr_raw = {}
for u, ps in by_iso.items():
    if "R" in ps:
        amr_raw[u] = "R"
    elif "S" in ps:
        amr_raw[u] = "S"
    else:
        amr_raw[u] = "UNKNOWN"
pheno = {}
with (RPT / "final_external_bdq_phenotype_validation_isolates.csv").open(
    newline="", encoding="utf-8-sig"
) as f:
    for r in csv.DictReader(f):
        pheno[r["unique_id"]] = pred(r.get("bdq_binary_phenotype"))
tbp_raw = {}
with open(CALLS, encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        v = (r.get("tbprofiler_bdq") or "").strip().upper()
        tbp_raw[r["id"]] = "R" if v == "R" else "UNKNOWN"
grade_map = {}
with (RPT / "who_bdq_per_row_benchmark.csv").open(
    newline="", encoding="utf-8-sig"
) as f:
    for r in csv.DictReader(f):
        grade_map[r["who_row_id"]] = r["who_truth"]
iso_cat = defaultdict(set)
with (RPT / "final_external_bdq_phenotype_validation_variant_calls.csv").open(
    newline="", encoding="utf-8-sig"
) as f:
    for r in csv.DictReader(f):
        wid = r["matched_who_row_id"]
        if wid and wid in grade_map and (grade_map[wid] in ("R", "S")):
            iso_cat[r["unique_id"]].add(grade_map[wid])
catalog = {u: "R" if "R" in iso_cat.get(u, set()) else "S" for u in pheno}
assert set(pheno.values()) <= {
    "R",
    "S",
}, "phenotype contains values outside R/S: %s" % sorted(
    set(pheno.values()) - {"R", "S"}
)
n_no_call = {
    "AMR": sum((1 for u in pheno if amr_raw.get(u) == "UNKNOWN")),
    "TBP": sum((1 for u in pheno if tbp_raw.get(u) != "R")),
}
amr = {u: "R" if amr_raw.get(u) == "R" else "S" for u in pheno}
tbp_full = {u: "R" if tbp_raw.get(u) == "R" else "S" for u in pheno}
methods_calls = {"CAT": catalog, "TBP": tbp_full, "AMR": amr}


def confusion(calls):
    tp = tn = fp = fn = 0
    for u, p in pheno.items():
        c = calls.get(u, "S")
        if p == "R":
            tp += c == "R"
            fn += c != "R"
        else:
            fp += c == "R"
            tn += c != "R"
    return (tp, tn, fp, fn)


def exact_mcnemar(a, b, sel):
    b_, c_ = (0, 0)
    for u, p in pheno.items():
        if p != sel:
            continue
        x, y = (a.get(u, "S"), b.get(u, "S"))
        if x == "R" and y != "R":
            b_ += 1
        elif x != "R" and y == "R":
            c_ += 1
    if b_ + c_ == 0:
        return (b_, c_, 0, 1.0)
    return (b_, c_, b_ + c_, float(binomtest(min(b_, c_), b_ + c_, 0.5).pvalue))


def holm(pvals):
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    adj = [1.0] * m
    prev = 0.0
    for rank, i in enumerate(order):
        a = pvals[i] * (m - rank)
        adj[i] = min(1.0, max(prev, a))
        prev = adj[i]
    return adj


def cochran_q(vectors, ids):
    k = len(vectors)
    col_tot = [sum((v.get(u, 0) for u in ids)) for _, v in vectors]
    row_tot = [sum((v.get(u, 0) for _, v in vectors)) for u in ids]
    num = (k - 1) * (k * sum((t * t for t in col_tot)) - sum(col_tot) ** 2)
    den = k * sum(row_tot) - sum((r * r for r in row_tot))
    q = num / den if den else 0.0
    pval = float(chi2.sf(q, k - 1))
    return (q, pval, col_tot)


conf = {k: confusion(v) for k, v in methods_calls.items()}
print("confusion (TP/TN/FP/FN):")
for k, v in conf.items():
    print(f"  {k}: {v}")
print(f"no-call counts: {n_no_call}")
pairs = [("AMR", "TBP"), ("AMR", "CAT"), ("TBP", "CAT")]
raw = {}
print("\nMcNemar exact (raw P; Holm-adjusted per endpoint):")
for sel, tag in (("R", "sensitivity"), ("S", "specificity")):
    pvals = []
    for a, b in pairs:
        bb, cc, nn, p = exact_mcnemar(methods_calls[a], methods_calls[b], sel)
        raw[a, b, tag] = p
        pvals.append(p)
    adj = holm(pvals)
    for (a, b), p, pa in zip(pairs, pvals, adj):
        raw[a, b, tag + "_adj"] = pa
        print(f"  {a} vs {b:3s} {tag:11s} rawP={p:.4g} holmP={pa:.4g}")
r_ids = [u for u in pheno if pheno[u] == "R"]
s_ids = [u for u in pheno if pheno[u] == "S"]
n_r, n_s = (len(r_ids), len(s_ids))
q_sens, qp_sens, _ = cochran_q(
    [
        (k, {u: 1 if v.get(u, "S") == "R" else 0 for u in r_ids})
        for k, v in methods_calls.items()
    ],
    r_ids,
)
q_spec, qp_spec, _ = cochran_q(
    [
        (k, {u: 1 if v.get(u, "S") == "R" else 0 for u in s_ids})
        for k, v in methods_calls.items()
    ],
    s_ids,
)
print(f"\nCochran Q sensitivity: Q={q_sens:.3f} df=2 P={qp_sens:.4g}")
print(f"Cochran Q specificity: Q={q_spec:.3f} df=2 P={qp_spec:.4g}")
methods = ["Direct WHO\ncatalogue matching", "TB-Profiler\nv6.7.0", "AMR-Hunter"]
C_CAT, C_TBP, C_AMR = ("#9aa5b1", "#4c72b0", "#d62728")
colors = [C_CAT, C_TBP, C_AMR]
order = ["CAT", "TBP", "AMR"]
tp = [conf[k][0] for k in order]
fp = [conf[k][2] for k in order]
tn = [conf[k][1] for k in order]
fn = [conf[k][3] for k in order]
sens = [t / n_r * 100 for t in tp]
fpr = [f / n_s * 100 for f in fp]
ppv = [t / (t + f) * 100 for t, f in zip(tp, fp)]


def wilson_ci(k, n, z=1.96):
    p = k / n
    denom = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (center - half, center + half)


sens_ci = [wilson_ci(t, n_r) for t in tp]
fpr_ci = [wilson_ci(f, n_s) for f in fp]
ppv_ci = [wilson_ci(t, t + f) for t, f in zip(tp, fp)]
p_sens = raw["AMR", "TBP", "sensitivity"]
p_sens_adj = raw["AMR", "TBP", "sensitivity_adj"]
p_spec = raw["AMR", "TBP", "specificity"]
p_spec_adj = raw["AMR", "TBP", "specificity_adj"]
cat_sens_vs_both = max(
    raw["AMR", "CAT", "sensitivity_adj"], raw["TBP", "CAT", "sensitivity_adj"]
)
cat_spec_vs_both = max(
    raw["AMR", "CAT", "specificity_adj"], raw["TBP", "CAT", "specificity_adj"]
)
plt.rcParams.update({"font.size": 9, "font.family": "DejaVu Sans"})
fig, axes = plt.subplots(1, 3, figsize=(9.8, 3.9))
x = range(3)


def bars(ax, vals, ylim, ylabel, labels, title, ci, fmt="%.1f%%"):
    ax.bar(x, vals, width=0.62, color=colors, edgecolor="white", linewidth=0.8)
    for i, (lo, hi) in enumerate(ci):
        ax.plot([i, i], [lo * 100, hi * 100], color="#333333", lw=1.1, zorder=5)
        ax.plot(
            [i - 0.09, i + 0.09],
            [lo * 100, lo * 100],
            color="#333333",
            lw=1.1,
            zorder=5,
        )
        ax.plot(
            [i - 0.09, i + 0.09],
            [hi * 100, hi * 100],
            color="#333333",
            lw=1.1,
            zorder=5,
        )
    for i, v in enumerate(vals):
        ax.text(
            i,
            ci[i][1] * 100 + 0.15,
            fmt % v,
            ha="center",
            va="bottom",
            fontsize=8,
            fontweight="bold",
            zorder=6,
        )
    ax.set_xticks(list(x))
    ax.set_xticklabels(methods, fontsize=8)
    ax.set_ylim(*ylim)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.set_title(title, fontsize=10, fontweight="bold", loc="left")
    ax.spines[["top", "right"]].set_visible(False)
    for i, lab in enumerate(labels):
        y_lab = max(0.12, (ylim[1] - ylim[0]) * 0.015)
        ax.text(
            i,
            y_lab,
            lab,
            ha="center",
            va="bottom",
            fontsize=6,
            color="white",
            alpha=0.9,
            zorder=6,
        )


def bracket(ax, x0, x1, y, text, tick_h=0.8):
    ax.plot([x0, x0, x1, x1], [y - tick_h, y, y, y - tick_h], color="black", lw=0.8)
    ax.text((x0 + x1) / 2, y + 0.2, text, ha="center", va="bottom", fontsize=7.5)


bars(
    axes[0],
    sens,
    (0, 57),
    "Sensitivity (%)",
    [f"{tp[i]}/{n_r}" for i in range(3)],
    "A",
    sens_ci,
)
if cat_sens_vs_both < 0.001:
    axes[0].text(0, 17.7, "***", ha="center", va="bottom", fontsize=10)
bracket(
    axes[0],
    1,
    2,
    55.2,
    (
        f"Holm-adjusted P = {p_sens_adj:.2f}"
        if p_sens_adj >= 0.05
        else f"* Holm-adjusted P = {p_sens_adj:.2f}"
    ),
)
bars(
    axes[1],
    fpr,
    (0, 2.6),
    "False-positive rate\n(FPR, %)",
    [f"{fp[i]}/{n_s}" for i in range(3)],
    "B",
    fpr_ci,
    fmt="%.2f%%",
)
if cat_spec_vs_both < 0.001:
    axes[1].text(0, 0.7, "***", ha="center", va="bottom", fontsize=10)
bracket(
    axes[1],
    1,
    2,
    2.1,
    (
        f"Holm-adjusted P = {p_spec_adj:.2f}"
        if p_spec_adj >= 0.05
        else f"* Holm-adjusted P = {p_spec_adj:.2f}"
    ),
    tick_h=0.35,
)
bars(
    axes[2],
    ppv,
    (0, 42),
    "Precision (PPV, %)",
    [f"{tp[i]}/{tp[i] + fp[i]}" for i in range(3)],
    "C",
    ppv_ci,
)
fig.tight_layout(w_pad=2.6)
for ext in ("png", "pdf", "svg"):
    fig.savefig(OUT_DIR / f"fig3_method_comparison.{ext}")
plt.close(fig)
print(f"\nFigure 3 (v3) written to {OUT_DIR}")
print("sens:", [round(v, 2) for v in sens])
print("fpr:", [round(v, 3) for v in fpr])
print("ppv:", [round(v, 2) for v in ppv])
print("TN:", tn, "FN:", fn)
