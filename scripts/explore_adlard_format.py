import csv
import io
import sys
import urllib.request
from collections import Counter, defaultdict

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
OUT = "d:\\python_work\\amr_hunter\\amr-hunter\\data\\external\\adlard2025"
RPT = "d:\\python_work\\amr_hunter\\amr-hunter\\reports\\publication"
who_path = OUT + "\\WHOv2_GARC1_RFUS.csv"
try:
    data = urllib.request.urlopen(
        "https://raw.githubusercontent.com/fowler-lab/tb-bdq-cat/main/catalogues/WHOv2/NC_000962.3_WHO-UCN-TB-2023.5_v2.0_GARC1_RFUS.csv",
        timeout=180,
    ).read()
    open(who_path, "wb").write(data)
    print("WHOv2 downloaded", len(data))
except Exception as e:
    print("WHOv2 dl err:", e)
    data = open(who_path, "rb").read()
rows = list(csv.reader(io.StringIO(data.decode("utf-8", "replace"))))
print("WHOv2 rows:", len(rows), "cols:", rows[0])
pred = Counter((r[8] for r in rows[1:]))
print("WHOv2 PREDICTION counts:", dict(pred))
for r in rows[1:12]:
    print("  WHOv2 row:", r[7], "->", r[8])
c2 = list(
    csv.reader(open(OUT + "\\catalogues_manuscript_catomatic_2.csv", encoding="utf-8"))
)
print(
    "\ncatomatic_2 rows:",
    len(c2),
    "pred counts:",
    dict(Counter((r[8] for r in c2[1:]))),
)
for r in c2[1:]:
    if r[8] == "R":
        print("  R row:", repr(r[7]), "evidence:", r[10][:80])
        if sum((1 for x in c2[1:] if x[8] == "R")) > 8:
            pass
vc = list(
    csv.DictReader(
        open(RPT + "\\virtual_mutation_candidates.csv", encoding="utf-8-sig")
    )
)
print("\nlocal candidates:", len(vc))
combo = Counter(
    ((r["candidate_role"], r["mutation_kind"], r["benchmark_role"]) for r in vc)
)
for k, v in sorted(combo.items()):
    print("  ", k, v)
