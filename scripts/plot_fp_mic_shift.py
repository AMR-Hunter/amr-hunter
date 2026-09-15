import csv
import os
import shutil
from collections import Counter, defaultdict
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

RPT = "d:\\python_work\\amr_hunter\\amr-hunter\\reports\\publication"
FIG = "d:\\python_work\\amr_hunter\\amr-hunter\\figures"
DESK = "C:\\Users\\zhang\\Desktop\\深度学习\\3、文章\\图"
UP = {"p.Ile67Leu", "p.Phe19Ser", "p.Ser68Asn"}
DN = {"p.Glu55Asp"}


def load():
    pheno = {}
    with open(
        RPT + "\\final_external_bdq_phenotype_validation_isolates.csv",
        encoding="utf-8-sig",
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
            wid = (r.get("matched_who_row_id") or "").strip()
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
            by[r["unique_id"]].append((g, m, vc, ev, p, wid))
    iso = {}
    for u, rs in by.items():
        ps = [x[4] for x in rs]
        iso[u] = "R" if "R" in ps else "S" if "S" in ps else "UNKNOWN"
    grade = {}
    with open(RPT + "\\who_bdq_per_row_benchmark.csv", encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            grade[r["who_row_id"]] = (r["who_grade"], r["who_truth"])
    return (pheno, by, iso, grade)


def micnum(s):
    if s == "0.5":
        return 0.5
    if s == "1.0":
        return 1.0
    if s == "2.0":
        return 2.0
    if s in (">2", ">1"):
        return 4.0
    return 0.0


def miccat(s):
    if s in (">2", ">1"):
        return 4.0
    try:
        return float(s.replace("<=", "").replace(">", ""))
    except ValueError:
        return 0.0


pheno, by, iso, grade = load()
tp = [u for u in pheno if pheno[u][0] == "R" and iso[u] == "R"]
fp = [u for u in pheno if pheno[u][0] == "S" and iso[u] == "R"]
fn = [u for u in pheno if pheno[u][0] == "R" and iso[u] != "R"]
assert len(tp) == 24 and len(fp) == 81 and (len(fn) == 47), (len(tp), len(fp), len(fn))
sens_rows = []
for name, predf in [
    (">=0.5 (standard)", lambda m: m >= 0.5),
    (">=1.0", lambda m: m >= 1.0),
    (">=2.0", lambda m: m >= 2.0),
]:
    tp_, fp_, fn_, tn_ = (0, 0, 0, 0)
    for u, (ph, mic) in pheno.items():
        rr = predf(micnum(mic))
        cc = iso[u] == "R"
        if rr:
            tp_ += cc
            fn_ += not cc
        else:
            tn_ += not cc
            fp_ += cc
    nR, nS = (tp_ + fn_, tn_ + fp_)
    sens_rows.append(
        [
            name,
            nR,
            nS,
            tp_,
            fp_,
            fn_,
            tn_,
            round(100 * tp_ / nR, 1),
            round(100 * tn_ / nS, 2),
            round(100 * tp_ / (tp_ + fp_), 1) if tp_ + fp_ else "",
            round(100 * tn_ / (tn_ + fn_), 2) if tn_ + fn_ else "",
        ]
    )
with open(
    RPT + "\\suppl_mic_cutoff_sensitivity.csv", "w", newline="", encoding="utf-8-sig"
) as f:
    wtr = csv.writer(f)
    wtr.writerow(
        [
            "resistance_definition",
            "n_resistant",
            "n_susceptible",
            "tp",
            "fp",
            "fn",
            "tn",
            "sensitivity_percent",
            "specificity_percent",
            "precision_percent",
            "npv_percent",
        ]
    )
    wtr.writerows(sens_rows)
print("=== MIC sensitivity table ===")
for r in sens_rows:
    print(r)
fnhi = [u for u in fn if micnum(pheno[u][1]) >= 1.0]
with open(
    RPT + "\\suppl_high_mic_fn_isolates.csv", "w", newline="", encoding="utf-8-sig"
) as f:
    wtr = csv.writer(f)
    wtr.writerow(["unique_id", "bdq_mic", "variants"])
    for u in sorted(fnhi, key=lambda x: -micnum(pheno[x][1])):
        vs = "; ".join(
            sorted(set(("%s:%s" % (g, m) for g, m, vc, ev, p, wid in by[u])))
        )
        wtr.writerow([u, pheno[u][1], vs])
print("=== high-MIC FN n =", len(fnhi), "===")
rows = []
for g, m in sorted({(x[0], x[1]) for u in fp + tp for x in by[u] if x[4] == "R"}):
    carriers = [u for u, rs in by.items() if any((x[0] == g and x[1] == m for x in rs))]
    n = len(carriers)
    nr = sum((1 for u in carriers if pheno[u][0] == "R"))
    whog = sorted(
        {
            grade[x[5]][0]
            for u in carriers
            for x in by[u]
            if x[0] == g
            and x[1] == m
            and (x[5] in grade)
            and (grade[x[5]][0] in ("G1", "G2"))
        }
    )
    frameshift = "fs" in m or m.endswith("*")
    who_cat = (
        ";".join(whog) if whog else "LoF-G1" if g == "Rv0678" and frameshift else "-"
    )
    mic = Counter((pheno[u][1] for u in carriers))
    rows.append(
        [
            g,
            m,
            n,
            nr,
            n - nr,
            round(100 * nr / n, 1),
            who_cat,
            "; ".join(("%s:%d" % (k, v) for k, v in sorted(mic.items()))),
        ]
    )
rows.sort(key=lambda x: -x[2])
with open(
    RPT + "\\suppl_variant_penetrance.csv", "w", newline="", encoding="utf-8-sig"
) as f:
    wtr = csv.writer(f)
    wtr.writerow(
        [
            "gene",
            "mutation",
            "carriers",
            "resistant",
            "susceptible",
            "PPV_percent",
            "WHO_grade",
            "MIC_distribution",
        ]
    )
    wtr.writerows(rows)
whom = set()
for u in fp:
    for g, m, vc, ev, p, wid in by[u]:
        if wid and grade.get(wid, ("",))[0] in ("G1", "G2"):
            whom.add(u)
        if g == "Rv0678" and ("fs" in m or m.endswith("*")):
            whom.add(u)
print("=== FP isolates with WHO-listed evidence:", len(whom), "of", len(fp), "===")
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
plt.rcParams.update(
    {"font.size": 10, "axes.spines.top": False, "axes.spines.right": False}
)
fig, axes = plt.subplots(
    1, 2, figsize=(9.4, 4.1), gridspec_kw={"width_ratios": [1.15, 1]}
)
ax = axes[0]
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
ax.set_title("A. False-positive MIC shift", loc="left", fontweight="bold")
ax.bar_label(b1, fmt="%.1f", fontsize=7.5, padding=1)
ax.bar_label(b2, fmt="%.1f", fontsize=7.5, padding=1)
ax.annotate(
    "42.0% vs 1.8%\nFisher P = 6.4 × 10$^{-36}$\nOR = 39.0",
    xy=(max(fpp[-1], allp[-1]) + 3, y[-1] + 0.19),
    fontsize=9,
    color="#8b0000",
    va="center",
    annotation_clip=False,
)
ax.set_xlim(0, max(fpp) * 1.55)
ax.legend(frameon=False, fontsize=8.5, loc="lower right")
ax.invert_yaxis()
ax = axes[1]
cuts = ["≥ 0.5", "≥ 1.0", "≥ 2.0"]
vals = [r[7] for r in sens_rows]
ns = ["24 / 71", "4 / 23", "1 / 11"]
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
ax.set_xlabel("MIC threshold defining resistance")
ax.set_title("B. Sensitivity vs MIC threshold", loc="left", fontweight="bold")
ax.annotate(
    "83.3% of true positives had MIC = 0.5 µg/mL",
    xy=(0.02, 0.94),
    xycoords="axes fraction",
    fontsize=8.5,
    color="#555555",
)
fig.tight_layout()
os.makedirs(FIG, exist_ok=True)
for ext in ("png", "pdf", "svg"):
    fig.savefig(FIG + f"\\cryptic_fp_mic_shift.{ext}", dpi=300, bbox_inches="tight")
    shutil.copy(
        FIG + f"\\cryptic_fp_mic_shift.{ext}", DESK + f"\\cryptic_fp_mic_shift.{ext}"
    )
print("figure saved")
