import csv
from collections import Counter, defaultdict

RPT = "d:\\python_work\\amr_hunter\\amr-hunter\\reports\\publication"
grade_map = {}
with open(RPT + "\\who_bdq_per_row_benchmark.csv", encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        grade_map[r["who_row_id"]] = (r["who_grade"], r["who_truth"])
iso_cat = defaultdict(set)
iso_pheno = {}
with open(
    RPT + "\\external_bdq_phenotype_validation_variant_calls.csv", encoding="utf-8-sig"
) as f:
    for r in csv.DictReader(f):
        uid = r["unique_id"]
        iso_pheno[uid] = r["bdq_binary_phenotype"]
        wid = r["matched_who_row_id"]
        if wid and wid in grade_map:
            g, truth = grade_map[wid]
            if truth == "R":
                iso_cat[uid].add("R")
            elif truth == "S":
                iso_cat[uid].add("S")
all_iso = {}
with open(
    RPT + "\\external_bdq_phenotype_validation_isolates.csv", encoding="utf-8-sig"
) as f:
    for r in csv.DictReader(f):
        all_iso[r["unique_id"]] = r["bdq_binary_phenotype"]


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
    return (tp, tn, fp, fn, unc, sen, spe, ppv)


calls1 = {}
for uid, s in iso_cat.items():
    if "R" in s:
        calls1[uid] = "R"
    elif "S" in s:
        calls1[uid] = "S"
    else:
        calls1[uid] = None
t1 = conf(calls1)
print("目录匹配(无条目=UNKNOWN): TP/TN/FP/FN/UNC=%d/%d/%d/%d/%d" % t1[:5])
print("  sens=%.3f spec=%.4f ppv=%.3f" % (t1[5], t1[6], t1[7]))
calls2 = {uid: "R" if "R" in s else "S" for uid, s in iso_cat.items()}
t2 = conf(calls2)
print("目录匹配(无条目=S): TP/TN/FP/FN=%d/%d/%d/%d" % t2[:4])
print("  sens=%.3f spec=%.4f ppv=%.3f" % (t2[5], t2[6], t2[7]))
fh = defaultdict(str)
with open(RPT + "\\cryptic_five_level_isolates.csv", encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        fh[r["unique_id"]] = r["five_level_label"]
callsR = {u: "R" if l in ("R", "Likely R") else "S" for u, l in fh.items()}
t3 = conf(callsR)
print("AMR-Hunter(R+Likely R): TP/TN/FP/FN=%d/%d/%d/%d" % t3[:4])
print("  sens=%.3f spec=%.4f ppv=%.3f" % (t3[5], t3[6], t3[7]))
callsHR = {u: "R" if l == "R" else "S" for u, l in fh.items()}
t4 = conf(callsHR)
print("AMR-Hunter(hard R): TP/TN/FP/FN=%d/%d/%d/%d" % t4[:4])
print("  sens=%.3f spec=%.4f ppv=%.3f" % (t4[5], t4[6], t4[7]))
print("phenotype R total:", sum((1 for p in all_iso.values() if p == "R")))
