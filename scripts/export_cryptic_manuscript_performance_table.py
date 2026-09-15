from __future__ import annotations
import argparse
import csv
from pathlib import Path

DEFAULT_BINARY_SUMMARY = Path(
    "reports/publication/external_bdq_phenotype_validation_summary.csv"
)
DEFAULT_FIVE_LEVEL_SUMMARY = Path(
    "reports/publication/cryptic_five_level_label_summary.csv"
)
DEFAULT_CV_TRADEOFF = Path("reports/publication/cryptic_supervised_cv_tradeoff.csv")
DEFAULT_SITE_HOLDOUT = Path(
    "reports/publication/cryptic_supervised_site_holdout_summary.csv"
)
DEFAULT_CSV_OUTPUT = Path(
    "reports/publication/cryptic_manuscript_performance_table.csv"
)
DEFAULT_MD_OUTPUT = Path("reports/publication/cryptic_manuscript_performance_table.md")
FIELDNAMES = [
    "table_section",
    "analysis",
    "evaluation_design",
    "call_policy",
    "model",
    "threshold_target",
    "total_isolates_or_predictions",
    "phenotype_r",
    "phenotype_s",
    "evaluated",
    "coverage",
    "unknown_no_call",
    "indeterminate",
    "tp",
    "tn",
    "fp",
    "fn",
    "recall",
    "specificity",
    "precision",
    "accuracy",
    "notes",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a manuscript-ready CRyPTIC BDQ performance comparison table."
    )
    parser.add_argument("--binary-summary", type=Path, default=DEFAULT_BINARY_SUMMARY)
    parser.add_argument(
        "--five-level-summary", type=Path, default=DEFAULT_FIVE_LEVEL_SUMMARY
    )
    parser.add_argument("--cv-tradeoff", type=Path, default=DEFAULT_CV_TRADEOFF)
    parser.add_argument("--site-holdout", type=Path, default=DEFAULT_SITE_HOLDOUT)
    parser.add_argument("--csv-output", type=Path, default=DEFAULT_CSV_OUTPUT)
    parser.add_argument("--md-output", type=Path, default=DEFAULT_MD_OUTPUT)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)


def as_int(value: str) -> int:
    return int(float(value)) if value else 0


def ratio(num: int, den: int) -> str:
    return f"{num / den:.4f}" if den else ""


def binary_row(summary_rows: list[dict[str, str]]) -> dict[str, str]:
    row = summary_rows[0]
    return {
        "table_section": "Frozen external validation",
        "analysis": "Frozen binary AMR-Hunter",
        "evaluation_design": "CRyPTIC HIGH-quality phenotype validation",
        "call_policy": "R/S isolate call; UNKNOWN counted in coverage denominator",
        "model": "frozen_binary_rules",
        "threshold_target": "",
        "total_isolates_or_predictions": row["total_isolates"],
        "phenotype_r": row["phenotype_r"],
        "phenotype_s": row["phenotype_s"],
        "evaluated": row["interpretable_isolates"],
        "coverage": row["coverage"],
        "unknown_no_call": row["unknown_isolates"],
        "indeterminate": "",
        "tp": row["tp"],
        "tn": row["tn"],
        "fp": row["fp"],
        "fn": row["fn"],
        "recall": row["sensitivity"],
        "specificity": row["specificity"],
        "precision": ratio(as_int(row["tp"]), as_int(row["tp"]) + as_int(row["fp"])),
        "accuracy": row["accuracy"],
        "notes": "Frozen result; accuracy excludes UNKNOWN/no-call isolates.",
    }


def five_level_row(
    label_rows: list[dict[str, str]],
    positive_labels: set[str],
    analysis: str,
    call_policy: str,
) -> dict[str, str]:
    total = sum((as_int(row["total_isolates"]) for row in label_rows))
    total_r = sum((as_int(row["phenotype_r"]) for row in label_rows))
    total_s = sum((as_int(row["phenotype_s"]) for row in label_rows))
    tp = sum(
        (
            as_int(row["phenotype_r"])
            for row in label_rows
            if row["label"] in positive_labels
        )
    )
    fp = sum(
        (
            as_int(row["phenotype_s"])
            for row in label_rows
            if row["label"] in positive_labels
        )
    )
    fn = total_r - tp
    tn = total_s - fp
    indeterminate = sum(
        (
            as_int(row["total_isolates"])
            for row in label_rows
            if row["label"] == "Indeterminate"
        )
    )
    return {
        "table_section": "Five-level interpretation layer",
        "analysis": analysis,
        "evaluation_design": "CRyPTIC phenotype validation with calibrated five-level labels",
        "call_policy": call_policy,
        "model": "five_level_rules",
        "threshold_target": "",
        "total_isolates_or_predictions": str(total),
        "phenotype_r": str(total_r),
        "phenotype_s": str(total_s),
        "evaluated": str(total),
        "coverage": "1.0000",
        "unknown_no_call": "0",
        "indeterminate": str(indeterminate),
        "tp": str(tp),
        "tn": str(tn),
        "fp": str(fp),
        "fn": str(fn),
        "recall": ratio(tp, tp + fn),
        "specificity": ratio(tn, tn + fp),
        "precision": ratio(tp, tp + fp),
        "accuracy": ratio(tp + tn, total),
        "notes": "Indeterminate isolates are retained in the table and treated as non-resistant calls for this binary projection.",
    }


def selected_supervised_row(
    rows: list[dict[str, str]],
    model: str,
    target: str,
    analysis: str,
    design: str,
    notes: str,
) -> dict[str, str]:
    matches = [
        row
        for row in rows
        if row["model"] == model and row["target_specificity"] == target
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one row for {model} at {target}, found {len(matches)}"
        )
    row = matches[0]
    tp = as_int(row["tp"])
    tn = as_int(row["tn"])
    fp = as_int(row["fp"])
    fn = as_int(row["fn"])
    total = as_int(row.get("total_predictions", row.get("total", "")))
    return {
        "table_section": "Supervised genotype-only upper bound",
        "analysis": analysis,
        "evaluation_design": design,
        "call_policy": "Training-fold threshold selected to satisfy target specificity",
        "model": model,
        "threshold_target": target,
        "total_isolates_or_predictions": str(total),
        "phenotype_r": str(tp + fn),
        "phenotype_s": str(tn + fp),
        "evaluated": str(total),
        "coverage": "1.0000",
        "unknown_no_call": "0",
        "indeterminate": "",
        "tp": str(tp),
        "tn": str(tn),
        "fp": str(fp),
        "fn": str(fn),
        "recall": row["recall"],
        "specificity": row["specificity"],
        "precision": row["precision"],
        "accuracy": row["accuracy"],
        "notes": notes,
    }


def build_rows(
    binary_summary: list[dict[str, str]],
    five_level_summary: list[dict[str, str]],
    cv_tradeoff: list[dict[str, str]],
    site_holdout: list[dict[str, str]],
) -> list[dict[str, str]]:
    return [
        binary_row(binary_summary),
        five_level_row(
            five_level_summary,
            {"R"},
            "Five-level hard R only",
            "Only hard R counted as resistant; Likely R/Indeterminate/Likely S/S counted as non-resistant",
        ),
        five_level_row(
            five_level_summary,
            {"R", "Likely R"},
            "Five-level R or Likely R",
            "R and Likely R counted as resistant; Indeterminate/Likely S/S counted as non-resistant",
        ),
        selected_supervised_row(
            cv_tradeoff,
            "no_secondary",
            "0.95",
            "Repeated 5x5 CV best high-specificity recall",
            "Repeated stratified 5-fold CV, aggregated across 5 repeats",
            "Upper-bound model; 42,675 predictions are 8,535 isolates x 5 repeats.",
        ),
        selected_supervised_row(
            cv_tradeoff,
            "no_mtrb",
            "0.90",
            "Repeated 5x5 CV relaxed-specificity recall ceiling",
            "Repeated stratified 5-fold CV, aggregated across 5 repeats",
            "Upper-bound model showing recall/specificity tradeoff.",
        ),
        selected_supervised_row(
            site_holdout,
            "no_secondary",
            "0.95",
            "Leave-one-site-out best high-specificity recall",
            "Leave-one-site-out validation using UNIQUEID site prefix",
            "Upper-bound model; target specificity is selected on training sites and may shift on held-out sites.",
        ),
        selected_supervised_row(
            site_holdout,
            "no_mtrb",
            "0.90",
            "Leave-one-site-out relaxed-specificity recall ceiling",
            "Leave-one-site-out validation using UNIQUEID site prefix",
            "Upper-bound model showing held-out site recall/specificity tradeoff.",
        ),
    ]


def markdown_table(rows: list[dict[str, str]]) -> str:
    columns = [
        "analysis",
        "evaluation_design",
        "model",
        "threshold_target",
        "coverage",
        "indeterminate",
        "tp",
        "fp",
        "fn",
        "recall",
        "specificity",
        "precision",
        "accuracy",
    ]
    lines = [
        "# CRyPTIC BDQ Manuscript Performance Table",
        "",
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join((row[column] for column in columns)) + " |")
    lines.extend(
        [
            "",
            "Notes:",
            "",
            "- Supervised rows are genotype-only upper-bound analyses, not frozen AMR-Hunter deployment claims.",
            "- Five-level rows are binary projections of calibrated labels; the indeterminate count remains explicit.",
            "- Repeated CV denominators count 8,535 isolates across 5 repeats.",
        ]
    )
    return "\n".join(lines) + "\n"


def write_markdown(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(markdown_table(rows), encoding="utf-8")


def main() -> int:
    args = parse_args()
    rows = build_rows(
        read_csv(args.binary_summary),
        read_csv(args.five_level_summary),
        read_csv(args.cv_tradeoff),
        read_csv(args.site_holdout),
    )
    write_csv(args.csv_output, rows)
    write_markdown(args.md_output, rows)
    print(f"csv={args.csv_output}")
    print(f"markdown={args.md_output}")
    print(f"rows={len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
