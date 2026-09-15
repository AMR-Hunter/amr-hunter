import csv
import io
import re
import sys
from collections import defaultdict

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
RPT = "d:\\python_work\\amr_hunter\\amr-hunter\\reports\\publication"
EXT = "d:\\python_work\\amr_hunter\\amr-hunter\\data\\external\\adlard2025"


def load_rows(path, mut_col, pred_col):
    with open(path, encoding="utf-8") as f:
        rows = list(csv.reader(f))[1:]
    out = []
    for r in rows:
        if len(r) > max(mut_col, pred_col):
            out.append((r[mut_col], r[pred_col]))
    return out


who = load_rows(EXT + "\\WHOv2_GARC1_RFUS.csv", 6, 7)
adl = load_rows(EXT + "\\catalogues_manuscript_catomatic_2.csv", 7, 8)
who_r = {m for m, p in who if p == "R"}
adl_r = {m for m, p in adl if p == "R"}
print("WHOv2 R rows:", len(who_r), "| Adlard 2a R rows:", len(adl_r))
print("WHOv2 sample:", sorted(who_r)[:6])
print("Adlard sample:", sorted(adl_r)[:6])
reg_who = [
    m
    for m in who_r
    if re.search("[-@]\\d+[A-Z]>", m) or re.search("[A-Z]\\d+[A-Z]$", m) is None
]
print(
    "WHOv2 non-aa rows:",
    [m for m in who_r if not re.fullmatch("[A-Za-z0-9_.]+@[A-Z]\\d+[A-Z!]", m)][:10],
)


def norm_who_like(gene, mut):
    m = re.fullmatch("([A-Z])(\\d+)([A-Z*])", mut or "")
    if m:
        aa = m.group(1) + m.group(2) + ("!" if m.group(3) == "*" else m.group(3))
        return f"{gene}@{aa}"
    m = re.fullmatch("RV0678_MMPL5_INTERGENIC:([A-Z])(\\d+)([A-Z])", mut or "")
    if m:
        return f"REG:{gene}@{m.group(2)}{m.group(1)}>{m.group(3)}"
    m = re.fullmatch("DEFAULT_AROUND_START:([A-Z])(\\d+)([A-Z])", mut or "")
    if m:
        return f"REG:{gene}@{m.group(2)}{m.group(1)}>{m.group(3)}"
    return f"{gene}@{mut}"


cand = list(
    csv.DictReader(
        open(RPT + "\\virtual_mutation_candidates.csv", encoding="utf-8-sig")
    )
)
res = [
    r
    for r in cand
    if r["resistance_label"] == "RESISTANT"
    and r["benchmark_role"] == "primary_bdq"
    and (r["gene"] in ("Rv0678", "atpE", "pepQ"))
]
print("\nlocal primary RESISTANT candidates:", len(res))
rows_out = []
confirmed_novel = []
for r in res:
    g = r["gene"]
    key = norm_who_like(g, r["mutation"])
    in_who = key in who_r
    in_adl = key in adl_r
    rows_out.append(
        [
            g,
            r["mutation"],
            r["mutation_kind"],
            r["final_status"],
            r["evidence_code"],
            key,
            "Y" if in_who else "",
            "Y" if in_adl else "",
        ]
    )
    if in_adl and (not in_who):
        confirmed_novel.append((g, r["mutation"], key))
with open(
    RPT + "\\adlard_candidate_comparison.csv", "w", newline="", encoding="utf-8-sig"
) as f:
    w = csv.writer(f)
    w.writerow(
        [
            "gene",
            "mutation",
            "mutation_kind",
            "final_status",
            "evidence_code",
            "normalized",
            "in_who2023_R",
            "in_adlard2025_R",
        ]
    )
    w.writerows(rows_out)
from collections import Counter

stat = Counter(((r[6], r[7]) for r in rows_out))
print("overlap stats (WHO, Adlard):", dict(stat))
print("\nAdlard-confirmed NOT in WHO (novel candidates independently supported):")
for g, m, k in confirmed_novel:
    print("  ", g, m, "->", k)
adl_stops = {m for m in adl_r if re.search("@[A-Z]\\d+!", m)}
print("\nAdlard stop-gain rows:", len(adl_stops))
