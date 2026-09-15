from __future__ import annotations
import argparse
from pathlib import Path
from export_who_bdq_lof_augmented_tables import (
    DEFAULT_CONFIG,
    DEFAULT_INPUT,
    DEFAULT_OUTPUT_DIR,
    augment_rows,
    build_gap_rows,
    build_stratified_rows,
    load_logic_types,
    read_rows,
    summarize_group,
    write_csv,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export unified WHO BDQ callable benchmark tables from frozen per-row CSV."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help="Frozen direct-match per-row WHO BDQ benchmark CSV.",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG,
        help="AMR-Hunter config.yaml with target gene logic_type values.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for unified publication CSV outputs.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    rows, input_fields = read_rows(args.input)
    logic_types = load_logic_types(args.config)
    unified = augment_rows(rows, logic_types)
    summary = [summarize_group("overall_primary", "all", unified)]
    stratified = build_stratified_rows(unified)
    gaps = build_gap_rows(unified)
    extra_fields = [
        "original_current_evaluator_scope",
        "original_amr_hunter_prediction",
        "original_matched_status",
        "original_correct",
        "original_gap_reason",
        "call_source",
        "consequence_inference",
        "inferred_functional_state",
        "inferred_logic_type",
    ]
    per_row_fields = input_fields + [
        field for field in extra_fields if field not in input_fields
    ]
    summary_fields = [
        "scope_type",
        "scope_value",
        "who_rows",
        "primary_evaluable_rows",
        "current_evaluator_scope_rows",
        "matched_rows",
        "direct_matched_rows",
        "lof_inferred_rows",
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
    write_csv(
        args.output_dir / "who_bdq_per_row_benchmark_unified.csv",
        unified,
        per_row_fields,
    )
    write_csv(
        args.output_dir / "who_bdq_benchmark_summary_unified.csv",
        summary,
        summary_fields,
    )
    write_csv(
        args.output_dir / "who_bdq_stratified_benchmark_unified.csv",
        stratified,
        summary_fields,
    )
    write_csv(args.output_dir / "who_bdq_gap_breakdown_unified.csv", gaps, gap_fields)
    overall = summary[0]
    print(f"per_row={args.output_dir / 'who_bdq_per_row_benchmark_unified.csv'}")
    print(f"summary={args.output_dir / 'who_bdq_benchmark_summary_unified.csv'}")
    print(
        "callable={matched_rows} direct={direct_matched_rows} lof_inferred={lof_inferred_rows} coverage={coverage} accuracy={matched_call_accuracy}".format(
            **overall
        )
    )
    print(f"remaining_unmatched={overall['unmatched_rows']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
