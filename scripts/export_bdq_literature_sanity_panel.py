from __future__ import annotations
import argparse
import csv
from pathlib import Path
from typing import Mapping, Sequence

DEFAULT_INPUT = Path("reports/publication/who_bdq_per_row_benchmark.csv")
DEFAULT_OUTPUT = Path("reports/publication/bdq_literature_positive_sanity_panel.csv")
WHO_CATALOGUE_URL = "https://www.who.int/publications/i/item/9789240082410"
ANDRIES_2014_URL = (
    "https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0102135"
)
PANEL_ROWS = [
    {
        "panel_role": "target_site_positive_control",
        "gene": "atpE",
        "mutation": "p.Ala63Pro",
        "mechanism_axis": "ATP synthase target-site mutation",
        "source_title": "WHO 2023 TB mutation catalogue; Andries et al. 2014 PLOS ONE",
        "source_url": f"{WHO_CATALOGUE_URL}; {ANDRIES_2014_URL}",
        "source_support": "WHO exact resistant catalogue row; literature supports atpE target-based BDQ resistance mechanisms.",
        "paper_use": "main text sanity-check example",
    },
    {
        "panel_role": "target_site_positive_control",
        "gene": "atpE",
        "mutation": "p.Ile66Met",
        "mechanism_axis": "ATP synthase target-site mutation",
        "source_title": "WHO 2023 TB mutation catalogue; Andries et al. 2014 PLOS ONE",
        "source_url": f"{WHO_CATALOGUE_URL}; {ANDRIES_2014_URL}",
        "source_support": "WHO exact resistant catalogue row; literature supports atpE target-based BDQ resistance mechanisms.",
        "paper_use": "supplementary sanity-check example",
    },
    {
        "panel_role": "target_site_positive_control",
        "gene": "atpE",
        "mutation": "p.Glu61Asp",
        "mechanism_axis": "ATP synthase target-site mutation",
        "source_title": "WHO 2023 TB mutation catalogue; Andries et al. 2014 PLOS ONE",
        "source_url": f"{WHO_CATALOGUE_URL}; {ANDRIES_2014_URL}",
        "source_support": "WHO exact resistant catalogue row; literature supports atpE target-based BDQ resistance mechanisms.",
        "paper_use": "supplementary sanity-check example",
    },
    {
        "panel_role": "rv0678_efflux_positive_control",
        "gene": "Rv0678",
        "mutation": "p.Gly121Arg",
        "mechanism_axis": "Rv0678/MmpS5-MmpL5 efflux regulation",
        "source_title": "WHO 2023 TB mutation catalogue; Andries et al. 2014 PLOS ONE",
        "source_url": f"{WHO_CATALOGUE_URL}; {ANDRIES_2014_URL}",
        "source_support": "WHO exact resistant catalogue row; literature links Rv0678 mutations to BDQ MIC increases and MmpS5-MmpL5 upregulation.",
        "paper_use": "main text sanity-check example",
    },
    {
        "panel_role": "rv0678_efflux_positive_control",
        "gene": "Rv0678",
        "mutation": "p.Leu117Arg",
        "mechanism_axis": "Rv0678/MmpS5-MmpL5 efflux regulation",
        "source_title": "WHO 2023 TB mutation catalogue; Andries et al. 2014 PLOS ONE",
        "source_url": f"{WHO_CATALOGUE_URL}; {ANDRIES_2014_URL}",
        "source_support": "WHO exact resistant catalogue row; literature links Rv0678 mutations to BDQ MIC increases and MmpS5-MmpL5 upregulation.",
        "paper_use": "supplementary sanity-check example",
    },
    {
        "panel_role": "rv0678_lof_positive_control",
        "gene": "Rv0678",
        "mutation": "p.Arg132*",
        "mechanism_axis": "Rv0678 loss-of-function / efflux derepression",
        "source_title": "WHO 2023 TB mutation catalogue; Andries et al. 2014 PLOS ONE",
        "source_url": f"{WHO_CATALOGUE_URL}; {ANDRIES_2014_URL}",
        "source_support": "WHO exact resistant catalogue row; literature supports Rv0678 disruption as BDQ resistance mechanism.",
        "paper_use": "main text LoF sanity-check example",
    },
    {
        "panel_role": "rv0678_gap_control",
        "gene": "Rv0678",
        "mutation": "p.Asp47fs",
        "mechanism_axis": "Rv0678 loss-of-function / efflux derepression",
        "source_title": "WHO 2023 TB mutation catalogue; Andries et al. 2014 PLOS ONE",
        "source_url": f"{WHO_CATALOGUE_URL}; {ANDRIES_2014_URL}",
        "source_support": "WHO exact resistant catalogue row; literature supports Rv0678 disruption as BDQ resistance mechanism.",
        "paper_use": "limitation example: known resistant frameshift not yet matched by final R/S export",
    },
    {
        "panel_role": "pepq_lof_positive_control",
        "gene": "pepQ",
        "mutation": "p.Gln5*",
        "mechanism_axis": "pepQ loss-of-function catalogue association",
        "source_title": "WHO 2023 TB mutation catalogue",
        "source_url": WHO_CATALOGUE_URL,
        "source_support": "WHO exact resistant catalogue row; use as catalogue-positive sanity check rather than independent phenotype validation.",
        "paper_use": "supplementary sanity-check example",
    },
    {
        "panel_role": "pepq_gap_control",
        "gene": "pepQ",
        "mutation": "LoF",
        "mechanism_axis": "pepQ generic loss-of-function catalogue category",
        "source_title": "WHO 2023 TB mutation catalogue",
        "source_url": WHO_CATALOGUE_URL,
        "source_support": "WHO exact resistant catalogue row; generic LoF category remains a documented current export gap.",
        "paper_use": "limitation example: generic LoF category not directly executable as a single variant",
    },
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export a small BDQ literature/WHO-positive sanity-check panel."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def read_who_rows(path: Path) -> dict[tuple[str, str], Mapping[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return {(row["gene"], row["mutation"]): row for row in rows}


def write_csv(
    path: Path, rows: Sequence[Mapping[str, object]], fields: Sequence[str]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def build_panel(
    who_rows: Mapping[tuple[str, str], Mapping[str, str]],
) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    missing: list[str] = []
    for item in PANEL_ROWS:
        key = (item["gene"], item["mutation"])
        who = who_rows.get(key)
        if who is None:
            missing.append(f"{key[0]} {key[1]}")
            continue
        output.append(
            {
                **item,
                "who_row_id": who["who_row_id"],
                "who_grade": who["who_grade"],
                "who_truth": who["who_truth"],
                "mutation_class": who["mutation_class"],
                "benchmark_role": who["benchmark_role"],
                "amr_hunter_prediction": who["amr_hunter_prediction"],
                "matched_status": who["matched_status"],
                "correct": who["correct"],
                "gap_reason": who["gap_reason"],
                "evidence_code": who["evidence_code"],
                "source_statuses": who["source_statuses"],
            }
        )
    if missing:
        raise RuntimeError(
            "Missing curated panel rows in WHO table: " + ", ".join(missing)
        )
    return output


def main() -> int:
    args = parse_args()
    who_rows = read_who_rows(args.input)
    panel = build_panel(who_rows)
    fields = [
        "panel_role",
        "gene",
        "mutation",
        "mechanism_axis",
        "who_row_id",
        "who_grade",
        "who_truth",
        "mutation_class",
        "benchmark_role",
        "amr_hunter_prediction",
        "matched_status",
        "correct",
        "gap_reason",
        "evidence_code",
        "source_statuses",
        "source_title",
        "source_url",
        "source_support",
        "paper_use",
    ]
    write_csv(args.output, panel, fields)
    matched = sum((1 for row in panel if row["matched_status"] == "matched_final"))
    print(f"panel={args.output}")
    print(f"panel_rows={len(panel)}")
    print(f"matched_rows={matched}")
    print(f"gap_or_unmatched_rows={len(panel) - matched}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
