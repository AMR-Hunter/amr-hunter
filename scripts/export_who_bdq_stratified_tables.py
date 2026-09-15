from __future__ import annotations
import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Mapping, Sequence

DEFAULT_INPUT = Path("reports/publication/who_bdq_per_row_benchmark.csv")
DEFAULT_OUTPUT_DIR = Path("reports/publication")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export stratified WHO BDQ benchmark and gap tables."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help="Frozen per-row WHO BDQ benchmark CSV.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for stratified publication CSV outputs.",
    )
    return parser.parse_args()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(
    path: Path, rows: Sequence[Mapping[str, object]], fields: Sequence[str]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def is_primary_evaluable(row: Mapping[str, str]) -> bool:
    return row["who_truth"] in {"R", "S"} and row["benchmark_role"] == "primary"


def is_current_scope(row: Mapping[str, str]) -> bool:
    return is_primary_evaluable(row) and row["current_evaluator_scope"] == "true"


def is_matched(row: Mapping[str, str]) -> bool:
    return is_primary_evaluable(row) and row["matched_status"] == "matched_final"


def compact_counter(counter: Counter[str], limit: int = 8) -> str:
    return ";".join((f"{key}:{value}" for key, value in counter.most_common(limit)))


def summarize_group(
    scope_type: str, scope_value: str, rows: Sequence[Mapping[str, str]]
) -> dict[str, object]:
    primary = [row for row in rows if is_primary_evaluable(row)]
    current_scope = [row for row in rows if is_current_scope(row)]
    matched = [row for row in rows if is_matched(row)]
    correct = [row for row in matched if row["correct"] == "true"]
    wrong = [row for row in matched if row["correct"] == "false"]
    unmatched = [row for row in primary if row["matched_status"] == "unmatched"]
    cm = Counter(((row["who_truth"], row["amr_hunter_prediction"]) for row in matched))
    gap_counts = Counter((row["gap_reason"] for row in unmatched if row["gap_reason"]))
    return {
        "scope_type": scope_type,
        "scope_value": scope_value,
        "who_rows": len(rows),
        "primary_evaluable_rows": len(primary),
        "current_evaluator_scope_rows": len(current_scope),
        "matched_rows": len(matched),
        "correct_rows": len(correct),
        "wrong_rows": len(wrong),
        "unmatched_rows": len(unmatched),
        "coverage": f"{len(matched) / len(primary):.4f}" if primary else "",
        "current_scope_coverage": (
            f"{len(matched) / len(current_scope):.4f}" if current_scope else ""
        ),
        "matched_call_accuracy": (
            f"{len(correct) / len(matched):.4f}" if matched else ""
        ),
        "tp": cm["R", "R"],
        "tn": cm["S", "S"],
        "fp": cm["S", "R"],
        "fn": cm["R", "S"],
        "gap_reasons": compact_counter(gap_counts),
    }


def group_by(
    rows: Iterable[Mapping[str, str]], keys: Sequence[str]
) -> dict[tuple[str, ...], list[Mapping[str, str]]]:
    groups: dict[tuple[str, ...], list[Mapping[str, str]]] = defaultdict(list)
    for row in rows:
        groups[tuple((row[key] or "BLANK" for key in keys))].append(row)
    return groups


def build_stratified_rows(rows: Sequence[Mapping[str, str]]) -> list[dict[str, object]]:
    output: list[dict[str, object]] = [summarize_group("overall_primary", "all", rows)]
    for key, items in sorted(group_by(rows, ["gene"]).items()):
        output.append(summarize_group("gene", key[0], items))
    primary_rows = [row for row in rows if row["benchmark_role"] == "primary"]
    for key, items in sorted(group_by(primary_rows, ["mutation_class"]).items()):
        output.append(summarize_group("mutation_class", key[0], items))
    for key, items in sorted(
        group_by(primary_rows, ["gene", "mutation_class"]).items()
    ):
        output.append(
            summarize_group("gene_mutation_class", f"{key[0]}|{key[1]}", items)
        )
    return [
        row
        for row in output
        if int(row["who_rows"]) > 0 and int(row["primary_evaluable_rows"]) > 0
    ]


def build_gap_rows(rows: Sequence[Mapping[str, str]]) -> list[dict[str, object]]:
    unmatched = [
        row
        for row in rows
        if is_primary_evaluable(row) and row["matched_status"] != "matched_final"
    ]
    groups = group_by(
        unmatched, ["gap_reason", "gene", "mutation_class", "who_truth", "who_grade"]
    )
    output: list[dict[str, object]] = []
    for key, items in sorted(groups.items()):
        gap_reason, gene, mutation_class, who_truth, who_grade = key
        examples = list(items[:5])
        output.append(
            {
                "gap_reason": gap_reason,
                "gene": gene,
                "mutation_class": mutation_class,
                "who_truth": who_truth,
                "who_grade": who_grade,
                "rows": len(items),
                "example_who_row_ids": "|".join(
                    (row["who_row_id"] for row in examples)
                ),
                "example_mutations": "|".join((row["mutation"] for row in examples)),
            }
        )
    return output


def main() -> int:
    args = parse_args()
    rows = read_rows(args.input)
    stratified = build_stratified_rows(rows)
    gaps = build_gap_rows(rows)
    stratified_fields = [
        "scope_type",
        "scope_value",
        "who_rows",
        "primary_evaluable_rows",
        "current_evaluator_scope_rows",
        "matched_rows",
        "correct_rows",
        "wrong_rows",
        "unmatched_rows",
        "coverage",
        "current_scope_coverage",
        "matched_call_accuracy",
        "tp",
        "tn",
        "fp",
        "fn",
        "gap_reasons",
    ]
    gap_fields = [
        "gap_reason",
        "gene",
        "mutation_class",
        "who_truth",
        "who_grade",
        "rows",
        "example_who_row_ids",
        "example_mutations",
    ]
    stratified_path = args.output_dir / "who_bdq_stratified_benchmark.csv"
    gap_path = args.output_dir / "who_bdq_gap_breakdown.csv"
    write_csv(stratified_path, stratified, stratified_fields)
    write_csv(gap_path, gaps, gap_fields)
    print(f"stratified={stratified_path}")
    print(f"gaps={gap_path}")
    print(f"stratified_rows={len(stratified)}")
    print(f"gap_rows={len(gaps)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
