import csv
import io
import re
import sys
import urllib.request

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
BASE = "https://raw.githubusercontent.com/fowler-lab/tb-bdq-cat/main/"
OUT = "d:\\python_work\\amr_hunter\\amr-hunter\\data\\external\\adlard2025"


def dl(path):
    data = urllib.request.urlopen(BASE + path, timeout=120).read()
    with open(OUT + "\\" + path.replace("/", "_"), "wb") as f:
        f.write(data)
    print("downloaded", path, len(data), "bytes")
    return data


readme = dl("README.md").decode("utf-8", "replace")
print("=== README ===")
print(readme[:2500])
xml = open(OUT + "\\fulltext.xml", encoding="utf-8").read()
txt = re.sub("<[^>]+>", " ", xml)
i = txt.find("catomatic")
while i >= 0 and i < len(txt):
    print("...", " ".join(txt[max(0, i - 150) : i + 200].split()))
    i = txt.find("catomatic", i + 1)
    if i > 40000:
        break
for p in [
    "catalogues/manuscript/catomatic_2.csv",
    "catalogues/manuscript/catomatic_3.csv",
    "catalogues/manuscript/catomatic_4.csv",
    "catalogues/manuscript/catomatic_5.csv",
    "catalogues/manuscript/catomatic_6.csv",
    "catalogues/manuscript/catomatic_1.csv",
]:
    try:
        data = dl(p)
    except Exception as e:
        print("ERR", p, e)
        continue
    rows = list(csv.reader(io.StringIO(data.decode("utf-8", "replace"))))
    print("===", p, "rows:", len(rows), "cols:", rows[0])
    for r in rows[1:4]:
        print("   ", r[:12])
