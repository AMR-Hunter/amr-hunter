import csv
import re
import sqlite3
import sys
from collections import Counter

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
DB = "C:\\Users\\zhang\\Desktop\\深度学习\\3、文章\\内部数据\\outdir\\analysis\\amr_hunter_db_final_20260908.db"
RPT = "d:\\python_work\\amr_hunter\\amr-hunter\\reports\\publication"
OUT = RPT + "\\saturation_R_tiers_20260914.csv"
GID = {3: "Rv0678", 5: "atpE", 7: "pepQ"}
LOF_GENES = {"Rv0678", "pepQ"}
THREE2ONE = dict(
    zip(
        [
            "Ala",
            "Arg",
            "Asn",
            "Asp",
            "Cys",
            "Gln",
            "Glu",
            "Gly",
            "His",
            "Ile",
            "Leu",
            "Lys",
            "Met",
            "Phe",
            "Pro",
            "Ser",
            "Thr",
            "Trp",
            "Tyr",
            "Val",
        ],
        "ARNDCQEGHILKMFPSTWYV",
    )
)


def norm_p(mutation: str):
    m = (mutation or "").strip()
    if not m.startswith("p."):
        return None
    m = m[2:]
    if "fs" in m:
        return None
    mm = re.fullmatch("([A-Z][a-z]{2})(\\d+)([A-Z][a-z]{2}|\\*|Ter|Xaa)", m)
    if not mm:
        return None
    a, b, d = mm.groups()
    if a not in THREE2ONE:
        return None
    d1 = "*" if d in ("*", "Ter") else THREE2ONE.get(d)
    if d1 is None:
        return None
    return THREE2ONE[a] + b + d1


def change_class(aa: str) -> str:
    m = re.fullmatch("([A-Z*])(\\d+)([A-Z*])", aa or "")
    if not m:
        return "other"
    a, b, d = m.groups()
    if a == "*":
        return "stop_lost"
    if d == "*":
        return "stop_gained"
    return "missense"


def main() -> int:
    con = sqlite3.connect(DB)
    ours = set()
    for gid, pos, ref, alt, aa in con.execute(
        "SELECT gene_id, pos, ref, alt, aa_change FROM mutations\n           WHERE resistance_phenotype='RESISTANT' AND region_type='CDS'\n             AND gene_id IN (3,5,7)"
    ):
        ours.add((GID[gid], pos, ref, alt, (aa or "").strip()))
    con.close()
    rows = list(
        csv.DictReader(
            open(RPT + "\\who_bdq_per_row_benchmark.csv", encoding="utf-8-sig")
        )
    )
    who_R, who_G3 = (set(), set())
    row_kind = Counter()
    for r in rows:
        mut, eff = (r["mutation"], r["effect"])
        if r["who_truth"] == "R":
            if mut == "LoF" or eff == "LoF":
                row_kind["generic_LoF"] += 1
            elif "fs" in mut or eff == "frameshift":
                row_kind["frameshift"] += 1
            elif mut.endswith("*") or eff == "stop_gained":
                row_kind["stop_gained"] += 1
            elif re.match("p\\.(Met1|M1)\\?", mut):
                row_kind["initiation_codon"] += 1
            else:
                row_kind["missense"] += 1
        aa = norm_p(mut)
        if not aa:
            continue
        if r["who_truth"] == "R":
            who_R.add((r["gene"], aa))
        elif r["who_truth"] == "":
            who_G3.add((r["gene"], aa))
    tier1a = {v for v in ours if (v[0], v[4]) in who_R}
    tier1b = {
        v
        for v in ours
        if v not in tier1a
        and change_class(v[4]) == "stop_gained"
        and (v[0] in LOF_GENES)
    }
    tier2 = {
        v
        for v in ours
        if v not in tier1a and v not in tier1b and ((v[0], v[4]) in who_G3)
    }
    tier3 = ours - tier1a - tier1b - tier2
    matched_rows = len({(v[0], v[4]) for v in tier1a})
    n_t1 = len(tier1a) + len(tier1b)
    assert len(ours) == 1643, len(ours)
    assert n_t1 + len(tier2) + len(tier3) == len(ours)
    print(f"universe (CDS RESISTANT, 10 target genes)          : {len(ours)}")
    print(
        f"  Tier 1a exact grade 1-2 match             : {len(tier1a)}  ({matched_rows}/{sum((v for k, v in row_kind.items() if k in ('missense', 'stop_gained')))} SNV-resolvable catalogue R rows matched)"
    )
    print(f"  Tier 1b nonsense covered by generic LoF   : {len(tier1b)}")
    print(f"  Tier 1 total (catalogue-anchored)         : {n_t1}")
    print(f"  Tier 2 grade-3 listed only                : {len(tier2)}")
    print(f"  Tier 3 absent from catalogue              : {len(tier3)}")
    print(
        f"  Tier 3 composition                        : {dict(Counter((change_class(v[4]) for v in tier3)))}"
    )
    print(f"\nWHO R rows by kind: {dict(row_kind)}  (total {sum(row_kind.values())})")
    lof_inferred = (
        row_kind["frameshift"] + row_kind["generic_LoF"] + row_kind["initiation_codon"]
    )
    print(
        f"LoF-inferred entries in the 461-entry concordance analysis: {lof_inferred} (frameshift + generic LoF + initiation codon)"
    )
    with open(OUT, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["item", "value", "definition"])
        w.writerow(
            [
                "universe",
                len(ours),
                "final DB, RESISTANT + CDS, 10 target genes, distinct (gene,pos,ref,alt)",
            ]
        )
        w.writerow(
            [
                "tier1a_exact_grade12_match",
                len(tier1a),
                f"aa change equals a WHO grade 1-2 entry (p. notation normalised); matches {matched_rows} catalogue rows",
            ]
        )
        w.writerow(
            [
                "tier1b_nonsense_generic_LoF",
                len(tier1b),
                "stop-gained variant in Rv0678/pepQ covered by the catalogue's generic LoF entries",
            ]
        )
        w.writerow(["tier1_catalogue_anchored", n_t1, "tier1a + tier1b"])
        w.writerow(
            ["tier2_grade3_listed", len(tier2), "matches only a WHO grade-3 entry"]
        )
        w.writerow(
            [
                "tier3_catalogue_absent",
                len(tier3),
                "no catalogue row match and not class-covered",
            ]
        )
        for k, v in row_kind.items():
            w.writerow([f"who_R_rows_{k}", v, "catalogue accounting, 88 R rows total"])
        w.writerow(["who_R_rows_total", sum(row_kind.values()), ""])
        w.writerow(
            [
                "lof_inferred_concordance_entries",
                lof_inferred,
                "frameshift + generic LoF + initiation codon rows (the 58 LoF-inferred entries)",
            ]
        )
    print(f"\nsaved: {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
