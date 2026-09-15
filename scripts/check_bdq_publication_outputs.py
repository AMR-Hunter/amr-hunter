from __future__ import annotations
import csv
import sys
from pathlib import Path
from typing import Iterable, Mapping, Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
from scripts.export_virtual_mutation_publication_tables import (
    candidate_rows,
    read_csv as read_virtual_csv,
    regulatory_assessment_rows,
    summarize,
)
from scripts.export_who_bdq_stratified_tables import (
    build_gap_rows,
    build_stratified_rows,
)
from scripts.export_bdq_literature_sanity_panel import build_panel
from scripts.export_who_bdq_lof_augmented_tables import (
    augment_rows as augment_who_lof_rows,
    build_gap_rows as build_lof_gap_rows,
    build_stratified_rows as build_lof_stratified_rows,
    load_logic_types,
    summarize_group as summarize_lof_group,
)

EXPECTED_WHO = {
    "who_total_rows": 1465,
    "primary_evaluable": 441,
    "current_scope": 438,
    "matched": 382,
    "correct": 382,
    "wrong": 0,
    "unmatched": 59,
    "tp": 30,
    "tn": 352,
    "fp": 0,
    "fn": 0,
}
EXPECTED_WHO_GAPS = {
    "unsupported_non_snp_variant": 56,
    "generic_lof_category": 2,
    "regulatory_out_of_config_window": 1,
}
EXPECTED_WHO_STRATIFIED_ROWS = {
    "overall_primary": 1,
    "gene": 8,
    "mutation_class": 7,
    "gene_mutation_class": 21,
}
EXPECTED_VIRTUAL_ROWS = {
    "snapshot": 33482,
    "summary": 11,
    "candidates": 507,
    "regulatory": 10,
}
EXPECTED_LITERATURE_PANEL = {"rows": 9, "matched": 7, "unmatched": 2}
EXPECTED_WHO_LOF_AUGMENTED = {
    "primary_evaluable": 441,
    "current_scope": 440,
    "matched": 440,
    "direct_matched": 382,
    "lof_inferred": 58,
    "correct": 440,
    "wrong": 0,
    "unmatched": 1,
    "tp": 88,
    "tn": 352,
    "fp": 0,
    "fn": 0,
}
EXPECTED_WHO_LOF_AUGMENTED_GAPS = {"regulatory_out_of_config_window": 1}
EXTERNAL_VALIDATION_FILES = {
    "summary": "external_bdq_phenotype_validation_summary.csv",
    "isolates": "external_bdq_phenotype_validation_isolates.csv",
    "variant_calls": "external_bdq_phenotype_validation_variant_calls.csv",
}
EXTERNAL_SUMMARY_FIELDS = [
    "dataset",
    "filter_name",
    "total_isolates",
    "phenotype_r",
    "phenotype_s",
    "interpretable_isolates",
    "unknown_isolates",
    "coverage",
    "accuracy",
    "tp",
    "tn",
    "fp",
    "fn",
    "sensitivity",
    "specificity",
    "notes",
]
EXTERNAL_ISOLATE_FIELDS = [
    "dataset",
    "unique_id",
    "ena_sample",
    "bdq_binary_phenotype",
    "bdq_mic",
    "bdq_phenotype_quality",
    "amr_hunter_isolate_prediction",
    "correct",
    "call_source_summary",
    "top_resistance_variant",
    "interpretable_variant_count",
    "unknown_reason",
    "vcf_path",
]
EXTERNAL_VARIANT_FIELDS = [
    "dataset",
    "unique_id",
    "ena_sample",
    "gene",
    "mutation",
    "region_type",
    "variant_class",
    "amr_hunter_prediction",
    "functional_state",
    "evidence_code",
    "call_source",
    "matched_who_row_id",
    "bdq_binary_phenotype",
    "bdq_mic",
    "vcf_path",
]
EXTERNAL_ALLOWED_CALL_SOURCES = {
    "direct_record_match",
    "consequence_inferred",
    "snapshot_virtual",
    "unsupported_variant",
    "no_call",
}


def read_dict_csv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        return (rows, list(reader.fieldnames or []))


def normalize_rows(
    rows: Iterable[Mapping[str, object]], fieldnames: Sequence[str]
) -> list[dict[str, str]]:
    return [{field: str(row.get(field, "")) for field in fieldnames} for row in rows]


def check_equal(
    label: str, actual: object, expected: object, failures: list[str]
) -> None:
    if actual != expected:
        failures.append(f"{label}: expected {expected!r}, got {actual!r}")


def format_ratio(numerator: int, denominator: int) -> str:
    if denominator == 0:
        return ""
    return f"{numerator / denominator:.4f}"


def check_required_fields(
    label: str, fieldnames: Sequence[str], required: Sequence[str], failures: list[str]
) -> None:
    missing = [field for field in required if field not in fieldnames]
    if missing:
        failures.append(f"{label}: missing field(s) {', '.join(missing)}")


def check_who_outputs(publication_dir: Path, failures: list[str]) -> None:
    per_row, _ = read_dict_csv(publication_dir / "who_bdq_per_row_benchmark.csv")
    summary, _ = read_dict_csv(publication_dir / "who_bdq_benchmark_summary.csv")
    stratified_saved, stratified_fields = read_dict_csv(
        publication_dir / "who_bdq_stratified_benchmark.csv"
    )
    gaps_saved, gap_fields = read_dict_csv(
        publication_dir / "who_bdq_gap_breakdown.csv"
    )
    literature_saved, literature_fields = read_dict_csv(
        publication_dir / "bdq_literature_positive_sanity_panel.csv"
    )
    primary = [
        row
        for row in per_row
        if row["who_truth"] in {"R", "S"} and row["benchmark_role"] == "primary"
    ]
    current_scope = [row for row in primary if row["current_evaluator_scope"] == "true"]
    matched = [row for row in primary if row["matched_status"] == "matched_final"]
    correct = [row for row in matched if row["correct"] == "true"]
    wrong = [row for row in matched if row["correct"] == "false"]
    unmatched = [row for row in primary if row["matched_status"] == "unmatched"]
    metrics = {
        "who_total_rows": len(per_row),
        "primary_evaluable": len(primary),
        "current_scope": len(current_scope),
        "matched": len(matched),
        "correct": len(correct),
        "wrong": len(wrong),
        "unmatched": len(unmatched),
        "tp": sum(
            (
                1
                for row in matched
                if row["who_truth"] == "R" and row["amr_hunter_prediction"] == "R"
            )
        ),
        "tn": sum(
            (
                1
                for row in matched
                if row["who_truth"] == "S" and row["amr_hunter_prediction"] == "S"
            )
        ),
        "fp": sum(
            (
                1
                for row in matched
                if row["who_truth"] == "S" and row["amr_hunter_prediction"] == "R"
            )
        ),
        "fn": sum(
            (
                1
                for row in matched
                if row["who_truth"] == "R" and row["amr_hunter_prediction"] == "S"
            )
        ),
    }
    for key, expected in EXPECTED_WHO.items():
        check_equal(f"WHO {key}", metrics[key], expected, failures)
    gap_counts: dict[str, int] = {}
    for row in unmatched:
        gap = row["gap_reason"]
        gap_counts[gap] = gap_counts.get(gap, 0) + 1
    check_equal("WHO gap counts", gap_counts, EXPECTED_WHO_GAPS, failures)
    overall = next(
        (
            row
            for row in summary
            if row["scope_type"] == "overall_primary" and row["scope_value"] == "all"
        ),
        None,
    )
    if overall is None:
        failures.append("WHO summary missing overall_primary/all row")
        return
    check_equal(
        "WHO summary matched rows",
        overall["matched_rows"],
        str(EXPECTED_WHO["matched"]),
        failures,
    )
    check_equal(
        "WHO summary matched accuracy",
        overall["matched_call_accuracy"],
        "1.0000",
        failures,
    )
    check_equal("WHO summary coverage", overall["coverage"], "0.8662", failures)
    check_equal(
        "WHO summary current coverage",
        overall["current_scope_coverage"],
        "0.8721",
        failures,
    )
    stratified = normalize_rows(build_stratified_rows(per_row), stratified_fields)
    gaps = normalize_rows(build_gap_rows(per_row), gap_fields)
    if stratified != stratified_saved:
        failures.append("WHO stratified benchmark table differs from regenerated rows")
    if gaps != gaps_saved:
        failures.append("WHO gap breakdown table differs from regenerated rows")
    scope_counts: dict[str, int] = {}
    for row in stratified_saved:
        scope_type = row["scope_type"]
        scope_counts[scope_type] = scope_counts.get(scope_type, 0) + 1
    check_equal(
        "WHO stratified scope counts",
        scope_counts,
        EXPECTED_WHO_STRATIFIED_ROWS,
        failures,
    )
    check_equal("WHO gap breakdown rows", len(gaps_saved), 7, failures)
    literature = normalize_rows(
        build_panel({(row["gene"], row["mutation"]): row for row in per_row}),
        literature_fields,
    )
    if literature != literature_saved:
        failures.append(
            "BDQ literature-positive sanity panel differs from regenerated rows"
        )
    literature_matched = sum(
        (1 for row in literature_saved if row["matched_status"] == "matched_final")
    )
    check_equal(
        "BDQ literature panel rows",
        len(literature_saved),
        EXPECTED_LITERATURE_PANEL["rows"],
        failures,
    )
    check_equal(
        "BDQ literature panel matched rows",
        literature_matched,
        EXPECTED_LITERATURE_PANEL["matched"],
        failures,
    )
    check_equal(
        "BDQ literature panel unmatched rows",
        len(literature_saved) - literature_matched,
        EXPECTED_LITERATURE_PANEL["unmatched"],
        failures,
    )


def check_who_lof_augmented_outputs(
    root: Path, publication_dir: Path, failures: list[str]
) -> None:
    per_row, _ = read_dict_csv(publication_dir / "who_bdq_per_row_benchmark.csv")
    augmented_saved, augmented_fields = read_dict_csv(
        publication_dir / "who_bdq_per_row_benchmark_lof_augmented.csv"
    )
    summary_saved, summary_fields = read_dict_csv(
        publication_dir / "who_bdq_benchmark_summary_lof_augmented.csv"
    )
    stratified_saved, stratified_fields = read_dict_csv(
        publication_dir / "who_bdq_stratified_benchmark_lof_augmented.csv"
    )
    gaps_saved, gap_fields = read_dict_csv(
        publication_dir / "who_bdq_gap_breakdown_lof_augmented.csv"
    )
    logic_types = load_logic_types(root / "config" / "config.yaml")
    augmented = augment_who_lof_rows(per_row, logic_types)
    summary = [summarize_lof_group("overall_primary", "all", augmented)]
    stratified = build_lof_stratified_rows(augmented)
    gaps = build_lof_gap_rows(augmented)
    comparisons = [
        (
            "WHO LoF-augmented per-row table",
            normalize_rows(augmented, augmented_fields),
            augmented_saved,
        ),
        (
            "WHO LoF-augmented summary table",
            normalize_rows(summary, summary_fields),
            summary_saved,
        ),
        (
            "WHO LoF-augmented stratified table",
            normalize_rows(stratified, stratified_fields),
            stratified_saved,
        ),
        ("WHO LoF-augmented gap table", normalize_rows(gaps, gap_fields), gaps_saved),
    ]
    for label, actual, expected in comparisons:
        if actual != expected:
            failures.append(f"{label}: regenerated rows differ from committed CSV")
    overall = summary[0]
    metrics = {
        "primary_evaluable": int(overall["primary_evaluable_rows"]),
        "current_scope": int(overall["current_evaluator_scope_rows"]),
        "matched": int(overall["matched_rows"]),
        "direct_matched": int(overall["direct_matched_rows"]),
        "lof_inferred": int(overall["lof_inferred_rows"]),
        "correct": int(overall["correct_rows"]),
        "wrong": int(overall["wrong_rows"]),
        "unmatched": int(overall["unmatched_rows"]),
        "tp": int(overall["tp"]),
        "tn": int(overall["tn"]),
        "fp": int(overall["fp"]),
        "fn": int(overall["fn"]),
    }
    for key, expected in EXPECTED_WHO_LOF_AUGMENTED.items():
        check_equal(f"WHO LoF-augmented {key}", metrics[key], expected, failures)
    gap_counts: dict[str, int] = {}
    for row in gaps:
        gap_counts[str(row["gap_reason"])] = gap_counts.get(
            str(row["gap_reason"]), 0
        ) + int(row["rows"])
    check_equal(
        "WHO LoF-augmented gap counts",
        gap_counts,
        EXPECTED_WHO_LOF_AUGMENTED_GAPS,
        failures,
    )


def check_who_unified_outputs(
    root: Path, publication_dir: Path, failures: list[str]
) -> None:
    per_row, _ = read_dict_csv(publication_dir / "who_bdq_per_row_benchmark.csv")
    unified_saved, unified_fields = read_dict_csv(
        publication_dir / "who_bdq_per_row_benchmark_unified.csv"
    )
    summary_saved, summary_fields = read_dict_csv(
        publication_dir / "who_bdq_benchmark_summary_unified.csv"
    )
    stratified_saved, stratified_fields = read_dict_csv(
        publication_dir / "who_bdq_stratified_benchmark_unified.csv"
    )
    gaps_saved, gap_fields = read_dict_csv(
        publication_dir / "who_bdq_gap_breakdown_unified.csv"
    )
    logic_types = load_logic_types(root / "config" / "config.yaml")
    unified = augment_who_lof_rows(per_row, logic_types)
    summary = [summarize_lof_group("overall_primary", "all", unified)]
    stratified = build_lof_stratified_rows(unified)
    gaps = build_lof_gap_rows(unified)
    comparisons = [
        (
            "WHO unified per-row table",
            normalize_rows(unified, unified_fields),
            unified_saved,
        ),
        (
            "WHO unified summary table",
            normalize_rows(summary, summary_fields),
            summary_saved,
        ),
        (
            "WHO unified stratified table",
            normalize_rows(stratified, stratified_fields),
            stratified_saved,
        ),
        ("WHO unified gap table", normalize_rows(gaps, gap_fields), gaps_saved),
    ]
    for label, actual, expected in comparisons:
        if actual != expected:
            failures.append(f"{label}: regenerated rows differ from committed CSV")
    overall = summary[0]
    metrics = {
        "primary_evaluable": int(overall["primary_evaluable_rows"]),
        "current_scope": int(overall["current_evaluator_scope_rows"]),
        "matched": int(overall["matched_rows"]),
        "direct_matched": int(overall["direct_matched_rows"]),
        "lof_inferred": int(overall["lof_inferred_rows"]),
        "correct": int(overall["correct_rows"]),
        "wrong": int(overall["wrong_rows"]),
        "unmatched": int(overall["unmatched_rows"]),
        "tp": int(overall["tp"]),
        "tn": int(overall["tn"]),
        "fp": int(overall["fp"]),
        "fn": int(overall["fn"]),
    }
    for key, expected in EXPECTED_WHO_LOF_AUGMENTED.items():
        check_equal(f"WHO unified {key}", metrics[key], expected, failures)


def check_virtual_outputs(
    root: Path, publication_dir: Path, failures: list[str]
) -> None:
    snapshot = read_virtual_csv(root / "reports" / "amr_results_snapshot.csv")
    check_equal(
        "virtual snapshot rows",
        len(snapshot),
        EXPECTED_VIRTUAL_ROWS["snapshot"],
        failures,
    )
    summary = summarize(snapshot)
    candidates = candidate_rows(snapshot, limit_per_gene=25)
    regulatory = regulatory_assessment_rows(snapshot, root / "config" / "config.yaml")
    saved_summary, summary_fields = read_dict_csv(
        publication_dir / "virtual_mutation_summary.csv"
    )
    saved_candidates, candidate_fields = read_dict_csv(
        publication_dir / "virtual_mutation_candidates.csv"
    )
    saved_regulatory, regulatory_fields = read_dict_csv(
        publication_dir / "regulatory_example_assessment.csv"
    )
    check_equal(
        "virtual summary rows", len(summary), EXPECTED_VIRTUAL_ROWS["summary"], failures
    )
    check_equal(
        "virtual candidate rows",
        len(candidates),
        EXPECTED_VIRTUAL_ROWS["candidates"],
        failures,
    )
    check_equal(
        "virtual regulatory rows",
        len(regulatory),
        EXPECTED_VIRTUAL_ROWS["regulatory"],
        failures,
    )
    comparisons = [
        (
            "virtual summary table",
            normalize_rows(summary, summary_fields),
            saved_summary,
        ),
        (
            "virtual candidate table",
            normalize_rows(candidates, candidate_fields),
            saved_candidates,
        ),
        (
            "virtual regulatory table",
            normalize_rows(regulatory, regulatory_fields),
            saved_regulatory,
        ),
    ]
    for label, actual, expected in comparisons:
        if actual != expected:
            failures.append(f"{label}: regenerated rows differ from committed CSV")


def check_external_validation_outputs(
    publication_dir: Path, failures: list[str]
) -> bool:
    paths = {
        label: publication_dir / filename
        for label, filename in EXTERNAL_VALIDATION_FILES.items()
    }
    present = {label: path.exists() for label, path in paths.items()}
    if not any(present.values()):
        return False
    if not all(present.values()):
        missing = [label for label, exists in present.items() if not exists]
        failures.append(
            "External BDQ phenotype validation outputs are incomplete; missing "
            + ", ".join(missing)
        )
        return True
    summary_rows, summary_fields = read_dict_csv(paths["summary"])
    isolate_rows, isolate_fields = read_dict_csv(paths["isolates"])
    variant_rows, variant_fields = read_dict_csv(paths["variant_calls"])
    check_required_fields(
        "External validation summary table",
        summary_fields,
        EXTERNAL_SUMMARY_FIELDS,
        failures,
    )
    check_required_fields(
        "External validation isolate table",
        isolate_fields,
        EXTERNAL_ISOLATE_FIELDS,
        failures,
    )
    check_required_fields(
        "External validation variant table",
        variant_fields,
        EXTERNAL_VARIANT_FIELDS,
        failures,
    )
    if failures:
        return True
    if len(summary_rows) != 1:
        failures.append(
            f"External validation summary table: expected 1 row, got {len(summary_rows)}"
        )
        return True
    phenotype_labels = {row["bdq_binary_phenotype"] for row in isolate_rows}
    invalid_labels = sorted(
        (label for label in phenotype_labels if label not in {"R", "S"})
    )
    if invalid_labels:
        failures.append(
            "External validation isolate table has non-CRyPTIC R/S phenotype labels: "
            + ", ".join(invalid_labels)
        )
    interpretable = [
        row
        for row in isolate_rows
        if row["amr_hunter_isolate_prediction"] in {"R", "S"}
    ]
    unknown = [
        row for row in isolate_rows if row["amr_hunter_isolate_prediction"] == "UNKNOWN"
    ]
    missing_unknown_reason = [
        row["unique_id"] for row in unknown if not row["unknown_reason"]
    ]
    if missing_unknown_reason:
        failures.append(
            "External validation UNKNOWN isolates missing unknown_reason: "
            + ", ".join(missing_unknown_reason[:10])
        )
    total = len(isolate_rows)
    phenotype_r = sum((1 for row in isolate_rows if row["bdq_binary_phenotype"] == "R"))
    phenotype_s = sum((1 for row in isolate_rows if row["bdq_binary_phenotype"] == "S"))
    tp = sum(
        (
            1
            for row in interpretable
            if row["bdq_binary_phenotype"] == "R"
            and row["amr_hunter_isolate_prediction"] == "R"
        )
    )
    tn = sum(
        (
            1
            for row in interpretable
            if row["bdq_binary_phenotype"] == "S"
            and row["amr_hunter_isolate_prediction"] == "S"
        )
    )
    fp = sum(
        (
            1
            for row in interpretable
            if row["bdq_binary_phenotype"] == "S"
            and row["amr_hunter_isolate_prediction"] == "R"
        )
    )
    fn = sum(
        (
            1
            for row in interpretable
            if row["bdq_binary_phenotype"] == "R"
            and row["amr_hunter_isolate_prediction"] == "S"
        )
    )
    expected_summary = {
        "total_isolates": str(total),
        "phenotype_r": str(phenotype_r),
        "phenotype_s": str(phenotype_s),
        "interpretable_isolates": str(len(interpretable)),
        "unknown_isolates": str(len(unknown)),
        "coverage": format_ratio(len(interpretable), total),
        "accuracy": format_ratio(tp + tn, len(interpretable)),
        "tp": str(tp),
        "tn": str(tn),
        "fp": str(fp),
        "fn": str(fn),
        "sensitivity": format_ratio(tp, tp + fn),
        "specificity": format_ratio(tn, tn + fp),
    }
    summary = summary_rows[0]
    for key, expected in expected_summary.items():
        check_equal(
            f"External validation summary {key}", summary[key], expected, failures
        )
    for row in variant_rows:
        call_source = row["call_source"]
        if call_source not in EXTERNAL_ALLOWED_CALL_SOURCES:
            failures.append(
                f"External validation variant call has unsupported call_source {call_source!r}"
            )
        if call_source == "no_call" and row["amr_hunter_prediction"] != "UNKNOWN":
            failures.append(
                f"External validation no_call row {row['unique_id']} is not UNKNOWN"
            )
    return True


def main() -> int:
    publication_dir = REPO_ROOT / "reports" / "publication"
    failures: list[str] = []
    check_who_outputs(publication_dir, failures)
    check_who_lof_augmented_outputs(REPO_ROOT, publication_dir, failures)
    check_who_unified_outputs(REPO_ROOT, publication_dir, failures)
    check_virtual_outputs(REPO_ROOT, publication_dir, failures)
    external_present = check_external_validation_outputs(publication_dir, failures)
    if failures:
        print("BDQ publication output check failed:")
        for failure in failures:
            print(f"- {failure}")
        return 1
    print("BDQ publication output check passed.")
    print(
        "WHO unified callable calls: 440/441, matched-call accuracy: 100%, coverage: 99.8%"
    )
    print("WHO call sources: 382 direct database-backed; 58 LoF consequence-inferred")
    print(
        "WHO direct-match audit: 382/441, matched-call accuracy: 100%, coverage: 86.6%"
    )
    print("Virtual mutation snapshot rows: 33482; candidate rows: 507")
    if external_present:
        print("External CRyPTIC BDQ phenotype validation outputs checked")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
