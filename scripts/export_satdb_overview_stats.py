from __future__ import annotations
import argparse
import csv
import sqlite3
from pathlib import Path
from core.analyzer import is_synonymous_cds_variant

FINAL_GENES = (
    "atpE",
    "Rv0678",
    "lpqB",
    "mtrA",
    "mtrB",
    "mmpS5",
    "mmpL5",
    "Rv1979c",
    "pepQ",
    "glpK",
)


def export_stats(db_path: Path, output_csv: Path | None = None) -> list[dict]:
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    placeholders = ",".join(("?" for _ in FINAL_GENES))
    rows = conn.execute(
        f"\n        SELECT DISTINCT m.gene_id, m.region_type, m.region_name, m.pos,\n                        m.ref, m.alt, m.aa_change, g.name AS gene_name\n        FROM mutations m\n        JOIN genes g ON m.gene_id = g.id\n        WHERE g.name IN ({placeholders})\n        ",
        FINAL_GENES,
    ).fetchall()
    gene_meta = {
        row["name"]: (
            str(row["logic_type"] or ""),
            str(row["scenario"] or ""),
            len(row["sequence"] or ""),
        )
        for row in conn.execute(
            f"SELECT name, logic_type, scenario, sequence FROM genes WHERE name IN ({placeholders})",
            FINAL_GENES,
        )
    }
    total = len(rows)
    coding = sum((1 for r in rows if str(r["region_type"]).upper() == "CDS"))
    regulatory = total - coding
    synonymous = sum(
        (
            1
            for r in rows
            if str(r["region_type"]).upper() == "CDS"
            and is_synonymous_cds_variant(dict(r))
        )
    )
    per_gene = {}
    for r in rows:
        gene = r["gene_name"]
        entry = per_gene.setdefault(gene, {"total": 0, "coding": 0, "regulatory": 0})
        entry["total"] += 1
        if str(r["region_type"]).upper() == "CDS":
            entry["coding"] += 1
        else:
            entry["regulatory"] += 1
    summary = {
        "metric": "BDQ-SatDB overview",
        "total_variants": total,
        "coding_variants": coding,
        "coding_pct": round(100.0 * coding / total, 2) if total else None,
        "regulatory_variants": regulatory,
        "regulatory_pct": round(100.0 * regulatory / total, 2) if total else None,
        "synonymous_coding": synonymous,
        "synonymous_pct_of_coding": (
            round(100.0 * synonymous / coding, 2) if coding else None
        ),
    }
    print(f"total: {total}")
    print(f"coding: {coding} ({summary['coding_pct']}%)")
    print(f"regulatory: {regulatory} ({summary['regulatory_pct']}%)")
    print(
        f"synonymous: {synonymous} ({summary['synonymous_pct_of_coding']}% of coding)"
    )
    print()
    print("per-gene:")
    for gene in FINAL_GENES:
        e = per_gene.get(gene, {"total": 0, "coding": 0, "regulatory": 0})
        logic, scenario, cds_len = gene_meta.get(gene, ("", "", 0))
        print(
            f"  {gene}: total={e['total']} coding={e['coding']} regulatory={e['regulatory']} cds_len={cds_len} logic={logic} scenario={scenario}"
        )
    conn.close()
    if output_csv is not None:
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        with output_csv.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(
                [
                    "metric",
                    "total_variants",
                    "coding_variants",
                    "coding_pct",
                    "regulatory_variants",
                    "regulatory_pct",
                    "synonymous_coding",
                    "synonymous_pct_of_coding",
                ]
            )
            writer.writerow(
                [
                    summary["metric"],
                    summary["total_variants"],
                    summary["coding_variants"],
                    summary["coding_pct"],
                    summary["regulatory_variants"],
                    summary["regulatory_pct"],
                    summary["synonymous_coding"],
                    summary["synonymous_pct_of_coding"],
                ]
            )
            writer.writerow([])
            writer.writerow(
                [
                    "gene",
                    "total",
                    "coding",
                    "regulatory",
                    "cds_length_bp",
                    "logic_type",
                    "scenario",
                ]
            )
            for gene in FINAL_GENES:
                e = per_gene.get(gene, {"total": 0, "coding": 0, "regulatory": 0})
                logic, scenario, cds_len = gene_meta.get(gene, ("", "", 0))
                writer.writerow(
                    [
                        gene,
                        e["total"],
                        e["coding"],
                        e["regulatory"],
                        cds_len,
                        logic,
                        scenario,
                    ]
                )
        print(f"\nwrote {output_csv}")
    return [summary]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--db", type=Path, required=True, help="AMR-Hunter SQLite database path."
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("reports/publication/satdb_overview_stats.csv"),
        help="Optional output CSV path.",
    )
    args = parser.parse_args()
    export_stats(args.db, args.output)


if __name__ == "__main__":
    main()
