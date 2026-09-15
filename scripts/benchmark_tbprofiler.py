import csv
from collections import defaultdict

RPT = "d:\\python_work\\amr_hunter\\amr-hunter\\reports\\publication"
CALLS = "C:\\Users\\zhang\\Desktop\\bdq_calls.csv"
all_iso = {}
with open(
    RPT + "\\external_bdq_phenotype_validation_isolates.csv", encoding="utf-8-sig"
) as f:
    for r in csv.DictReader(f):
        all_iso[r["unique_id"]] = r["bdq_binary_phenotype"]
tbp = {}
with open(CALLS, encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        bdq = (r.get("tbprofiler_bdq") or "").strip()
        tbp[r["id"]] = bdq


def conf(calls):
    tp = tn = fp = fn = unc = 0
    for uid, pheno in all_iso.items():
        c = calls.get(uid)
        if c is None:
            unc += 1
            continue
        if pheno == "R":
            if c == "R":
                tp += 1
            else:
                fn += 1
        elif c == "R":
            fp += 1
        else:
            tn += 1
    sen = tp / (tp + fn) if tp + fn else 0
    spe = tn / (tn + fp) if tn + fp else 0
    ppv = tp / (tp + fp) if tp + fp else 0
    npv = tn / (tn + fn) if tn + fn else 0
    return (tp, tn, fp, fn, unc, sen, spe, ppv, npv)


calls_u = {u: "R" if b == "R" else "S" if b == "S" else None for u, b in tbp.items()}
t1 = conf(calls_u)
print("TB-Profiler(空=UNKNOWN): TP/TN/FP/FN/UNC=%d/%d/%d/%d/%d" % t1[:5])
print("  sens=%.4f spec=%.4f ppv=%.4f npv=%.4f" % (t1[5], t1[6], t1[7], t1[8]))
calls_s = {u: "R" if b == "R" else "S" for u, b in tbp.items()}
t2 = conf(calls_s)
print("TB-Profiler(空=S): TP/TN/FP/FN=%d/%d/%d/%d" % t2[:4])
print("  sens=%.4f spec=%.4f ppv=%.4f npv=%.4f" % (t2[5], t2[6], t2[7], t2[8]))
print()
print(
    "覆盖株数: %d / %d (%.1f%%)"
    % (len(tbp), len(all_iso), 100 * len(tbp) / len(all_iso))
)
print("phenotype R total:", sum((1 for p in all_iso.values() if p == "R")))
print("TB-Profiler 判 R 株数:", sum((1 for b in tbp.values() if b == "R")))
grade_map = {}
with open(RPT + "\\who_bdq_per_row_benchmark.csv", encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        grade_map[r["who_row_id"]] = (r["who_grade"], r["who_truth"])
iso_cat = defaultdict(set)
with open(
    RPT + "\\external_bdq_phenotype_validation_variant_calls.csv", encoding="utf-8-sig"
) as f:
    for r in csv.DictReader(f):
        uid = r["unique_id"]
        wid = r["matched_who_row_id"]
        if wid and wid in grade_map:
            g, truth = grade_map[wid]
            if truth == "R":
                iso_cat[uid].add("R")
            elif truth == "S":
                iso_cat[uid].add("S")
from scipy.stats import chi2


def mcn_pval(tab):
    b, c = (tab[0][1], tab[1][0])
    if b + c == 0:
        return 1.0
    stat = (b - c) ** 2 / (b + c)
    return float(chi2.sf(stat, 1))


sens_tab = [[0, 0], [0, 0]]
spec_tab = [[0, 0], [0, 0]]
for uid, pheno in all_iso.items():
    c = "R" if iso_cat.get(uid) and "R" in iso_cat[uid] else "S"
    t = tbp.get(uid)
    t = "R" if t == "R" else "S"
    if pheno == "R":
        sens_tab[0 if t == "R" else 1][0 if c == "R" else 1] += 1
    else:
        spec_tab[0 if t == "R" else 1][0 if c == "R" else 1] += 1
print()
print("灵敏度 McNemar (TB-Profiler vs 目录):", sens_tab, "P=%.4g" % mcn_pval(sens_tab))
print("特异度 McNemar (TB-Profiler vs 目录):", spec_tab, "P=%.4g" % mcn_pval(spec_tab))
fh = defaultdict(str)
with open(RPT + "\\cryptic_five_level_isolates.csv", encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        fh[r["unique_id"]] = r["five_level_label"]
sens_tab2 = [[0, 0], [0, 0]]
spec_tab2 = [[0, 0], [0, 0]]
for uid, pheno in all_iso.items():
    a = fh.get(uid)
    a = "R" if a in ("R", "Likely R") else "S"
    t = tbp.get(uid)
    t = "R" if t == "R" else "S"
    if pheno == "R":
        sens_tab2[0 if t == "R" else 1][0 if a == "R" else 1] += 1
    else:
        spec_tab2[0 if t == "R" else 1][0 if a == "R" else 1] += 1
print(
    "灵敏度 McNemar (TB-Profiler vs AMR-Hunter):",
    sens_tab2,
    "P=%.4g" % mcn_pval(sens_tab2),
)
print(
    "特异度 McNemar (TB-Profiler vs AMR-Hunter):",
    spec_tab2,
    "P=%.4g" % mcn_pval(spec_tab2),
)
