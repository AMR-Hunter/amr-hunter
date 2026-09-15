import csv
from collections import defaultdict

RPT = "d:\\python_work\\amr_hunter\\amr-hunter\\reports\\publication"
grade_map = {}
with open(RPT + "\\who_bdq_per_row_benchmark.csv", encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        grade_map[r["who_row_id"]] = r["who_truth"]
iso_cat = defaultdict(set)
with open(
    RPT + "\\external_bdq_phenotype_validation_variant_calls.csv", encoding="utf-8-sig"
) as f:
    for r in csv.DictReader(f):
        wid = r["matched_who_row_id"]
        if wid and wid in grade_map:
            iso_cat[r["unique_id"]].add(grade_map[wid])
catalog = {u: "R" if "R" in s else "S" for u, s in iso_cat.items()}
amr = {}
with open(RPT + "\\cryptic_five_level_isolates.csv", encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        amr[r["unique_id"]] = "R" if r["five_level_label"] in ("R", "Likely R") else "S"
pheno = {}
with open(
    RPT + "\\external_bdq_phenotype_validation_isolates.csv", encoding="utf-8-sig"
) as f:
    for r in csv.DictReader(f):
        pheno[r["unique_id"]] = r["bdq_binary_phenotype"]
b = c = 0
both = 0
for u, p in pheno.items():
    if u not in catalog or u not in amr:
        continue
    if p != "R":
        continue
    both += 1
    if catalog[u] == "R" and amr[u] != "R":
        b += 1
    elif catalog[u] != "R" and amr[u] == "R":
        c += 1
print("pheno-R paired: %d, catalog-only-R(b)=%d, amr-only-R(c)=%d" % (both, b, c))
from math import comb

n = b + c
k = min(b, c)
pval = 2 * sum((comb(n, i) for i in range(k + 1))) * 0.5**n
print("McNemar exact p (sens, paired pheno-R):", pval)
b = c = 0
both = 0
for u, p in pheno.items():
    if u not in catalog or u not in amr:
        continue
    if p == "R":
        continue
    both += 1
    if catalog[u] == "R" and amr[u] != "R":
        b += 1
    elif catalog[u] != "R" and amr[u] == "R":
        c += 1
n = b + c
k = min(b, c)
pval = 2 * sum((comb(n, i) for i in range(k + 1))) * 0.5**n
print("pheno-S paired: %d, catalog-only-R(b)=%d, amr-only-R(c)=%d" % (both, b, c))
print("McNemar exact p (spec, paired pheno-S):", pval)
