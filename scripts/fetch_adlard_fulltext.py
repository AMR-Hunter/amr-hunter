import re
import urllib.request

OUT = (
    "d:\\python_work\\amr_hunter\\amr-hunter\\data\\external\\adlard2025\\fulltext.xml"
)
url = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12177157/fullTextXML"
req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
data = urllib.request.urlopen(req, timeout=120).read().decode("utf-8", "replace")
with open(OUT, "w", encoding="utf-8") as f:
    f.write(data)
print("chars:", len(data))
hits = re.findall(
    "(?:github\\.com|figshare|zenodo|osf\\.io)[^\\s<\\\"'\\)]*", data, re.I
)
print("repo hits:", sorted(set(hits)))
i = data.lower().find("data availability")
if i >= 0:
    print("=== data availability ===")
    print(re.sub("<[^>]+>", " ", data[i : i + 2500]))
else:
    print("no data availability section found")
