import csv
import re
import sys
from collections import defaultdict

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
RPT = "d:\\python_work\\amr_hunter\\amr-hunter\\reports\\publication"
EXT = "d:\\python_work\\amr_hunter\\amr-hunter\\data\\external\\adlard2025"


def load(path, mc, pc):
    rows = list(csv.reader(open(path, encoding="utf-8")))[1:]
    return [(r[mc], r[pc]) for r in rows if len(r) > max(mc, pc)]


who = load(EXT + "\\WHOv2_GARC1_RFUS.csv", 6, 7)
adl = load(EXT + "\\catalogues_manuscript_catomatic_2.csv", 7, 8)
who_r = {m for m, p in who if p == "R"}
adl_r = {m for m, p in adl if p == "R"}


def in_who(gene, aa1):
    return f"{gene}@{aa1}" in who_r or (aa1.endswith("!") and f"{gene}@*!" in who_r)


def in_adl(gene, aa1):
    return f"{gene}@{aa1}" in adl_r


def normalize(gene, mut):
    m = re.fullmatch("([A-Z])(\\d+)([A-Z*])", mut or "")
    if m:
        return m.group(1) + m.group(2) + ("!" if m.group(3) == "*" else m.group(3))
    return None


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
print("candidates:", len(res))
rows, confirmed_novel, novel_both_missing = ([], [], [])
for r in res:
    g = r["gene"]
    kind = r["mutation_kind"]
    a = normalize(g, r["mutation"])
    if a is None:
        iw = ia = ""
        note = "regulatory (mapping deferred)"
    else:
        iw = "Y" if in_who(g, a) else ""
        ia = "Y" if in_adl(g, a) else ""
        note = ""
        if iw and ia:
            pass
        elif ia and (not iw):
            confirmed_novel.append((g, r["mutation"], f"{g}@{a}"))
        elif not iw and (not ia):
            novel_both_missing.append((g, r["mutation"], f"{g}@{a}"))
    rows.append(
        [
            g,
            r["mutation"],
            kind,
            r["final_status"],
            r["evidence_code"],
            f"{g}@{a}" if a else "REG",
            iw,
            ia,
            note,
        ]
    )
with open(
    RPT + "\\adlard_interim_comparison.csv", "w", newline="", encoding="utf-8-sig"
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
            "in_who2023",
            "in_adlard2025",
            "note",
        ]
    )
    w.writerows(rows)
nw = sum((1 for r in rows if r[6] == "Y"))
na = sum((1 for r in rows if r[7] == "Y"))
nb = sum((1 for r in rows if r[6] == "Y" and r[7] == "Y"))
print(f"candidates in WHO={nw}, in Adlard={na}, both={nb}")
print("Adlard-confirmed but NOT in WHO (n=%d):" % len(confirmed_novel))
for g, m, k in confirmed_novel:
    print("   ", g, m)
print("neither WHO nor Adlard (candidate novel, n=%d):" % len(novel_both_missing))
for g, m, k in novel_both_missing[:40]:
    print("   ", g, m)
UP = {"p.Ile67Leu", "p.Phe19Ser", "p.Ser68Asn"}
DN = {"p.Glu55Asp"}
aa3 = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "Stop": "*",
}


def hgvsp1(g, mut):
    m = re.fullmatch("p\\.([A-Z][a-z]{2})(\\d+)([A-Za-z]{2,3}|\\*|fs)", mut or "")
    if not m:
        return None
    alt = m.group(3)
    if alt == "*":
        a = "!"
    elif alt == "fs":
        a = "_fs"
    else:
        a = aa3.get(alt)
        if a is None:
            return None
    return f"{g}@{aa3[m.group(1)]}{m.group(2)}{a}"


def fs_nt_match(g, mut):
    m = re.fullmatch("p\\.([A-Z][a-z]{2})(\\d+)fs", mut or "")
    if not m:
        return None
    n = int(m.group(2))
    for pos in range(3 * (n - 1) + 1, 3 * n + 1):
        pats = {f"{g}@{pos}_ins_x", f"{g}@{pos}_del_x"}
        for a in adl_r:
            if re.match(f"{re.escape(g)}@{pos}_(ins|del|indel)", a):
                return a
    return None


pheno = {
    r["unique_id"]: (r["bdq_binary_phenotype"],)
    for r in csv.DictReader(
        open(
            RPT + "\\final_external_bdq_phenotype_validation_isolates.csv",
            encoding="utf-8-sig",
        )
    )
}
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
tp = [u for u in pheno if pheno[u][0] == "R" and iso[u] == "R"]
fpset = {(g, m) for u in fp for g, m, p in by[u] if p == "R"}
tpset = {(g, m) for u in tp for g, m, p in by[u] if p == "R"}
print("\nFP variants:", len(fpset), "| TP variants:", len(tpset))


def in_adl_any(g, m):
    a = hgvsp1(g, m)
    if a and a in adl_r:
        return a
    return fs_nt_match(g, m)


fp_adl = {k: in_adl_any(*k) for k in fpset if in_adl_any(*k)}
tp_adl = {k: in_adl_any(*k) for k in tpset if in_adl_any(*k)}
fp_adl_iso = sum((1 for u in fp if any((in_adl_any(g, m) for g, m, p in by[u]))))
tp_adl_iso = sum((1 for u in tp if any((in_adl_any(g, m) for g, m, p in by[u]))))
print("FP variants in Adlard:", len(fp_adl), "-> isolates:", fp_adl_iso, "/", len(fp))
print("TP variants in Adlard:", len(tp_adl), "-> isolates:", tp_adl_iso, "/", len(tp))
for (g, m), a in sorted(fp_adl.items()):
    print("   FP in Adlard:", g, m, "->", a)
for (g, m), a in sorted(tp_adl.items()):
    print("   TP in Adlard:", g, m, "->", a)
fn = [u for u in pheno if pheno[u][0] == "R" and iso[u] != "R"]
fn_adl = [u for u in fn if any((in_adl_any(g, m) for g, m, p in by[u]))]
print("\nFN isolates carrying Adlard-R variants:", len(fn_adl), "/", len(fn))
for u in fn_adl:
    hits = [(g, m, in_adl_any(g, m)) for g, m, p in by[u] if in_adl_any(g, m)]
    print("   ", u, pheno[u], hits)
