from __future__ import annotations
import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Mapping, Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
from core.consequence import resolve_consequence_level_decision

try:
    import yaml
except ImportError:
    yaml = None
DEFAULT_INPUT = Path("reports/publication/who_bdq_per_row_benchmark.csv")
DEFAULT_OUTPUT_DIR = Path("reports/publication")
DEFAULT_CONFIG = Path("config/config.yaml")
LOF_MUTATION_CLASSES = {"frameshift", "start_lost", "generic_lof", "nonsense"}
MATCHED_STATUSES = {"matched_final", "consequence_inferred"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export LoF-augmented WHO BDQ benchmark tables from frozen per-row CSV."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help="Frozen per-row WHO BDQ benchmark CSV.",
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
        help="Directory for augmented publication CSV outputs.",
    )
    parser.add_argument(
        "--ignore-rv0678-neutral-fs",
        action="store_true",
        help="Disable the Rv0678 neutral-frameshift calibration for a pure-rule catalogue comparison (no external phenotype calibration applied).",
    )
    return parser.parse_args()


def read_rows(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return (list(reader), list(reader.fieldnames or []))


def write_csv(
    path: Path, rows: Sequence[Mapping[str, object]], fields: Sequence[str]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def load_logic_types(config_path: Path) -> dict[str, str]:
    if yaml is None:
        return load_logic_types_without_yaml(config_path)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    logic_types: dict[str, str] = {}
    for item in config.get("target_genes") or []:
        gene = str(item.get("gene") or "").strip()
        logic_type = str(item.get("logic_type") or "").strip().upper()
        if gene and logic_type:
            logic_types[gene] = logic_type
    return logic_types


def load_logic_types_without_yaml(config_path: Path) -> dict[str, str]:
    logic_types: dict[str, str] = {}
    in_target_genes = False
    current_gene = ""
    current_logic = ""

    def flush() -> None:
        nonlocal current_gene, current_logic
        if current_gene and current_logic:
            logic_types[current_gene] = current_logic
        current_gene = ""
        current_logic = ""

    for raw_line in config_path.read_text(encoding="utf-8").splitlines():
        stripped = raw_line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if not raw_line.startswith((" ", "-")) and stripped.endswith(":"):
            if in_target_genes:
                flush()
            in_target_genes = stripped == "target_genes:"
            continue
        if not in_target_genes:
            continue
        if stripped.startswith("- "):
            flush()
            stripped = stripped[2:].strip()
        if stripped.startswith("gene:"):
            current_gene = stripped.split(":", 1)[1].strip().strip("'\"")
        elif stripped.startswith("logic_type:"):
            current_logic = stripped.split(":", 1)[1].strip().strip("'\"").upper()
    if in_target_genes:
        flush()
    return logic_types


def can_infer_lof(row: Mapping[str, str], logic_types: Mapping[str, str]) -> bool:
    if row.get("who_truth") not in {"R", "S"}:
        return False
    if row.get("benchmark_role") != "primary":
        return False
    if row.get("matched_status") == "matched_final":
        return False
    if row.get("mutation_class") not in LOF_MUTATION_CLASSES:
        return False
    return bool(logic_types.get(str(row.get("gene") or "")))


def augment_rows(
    rows: Sequence[Mapping[str, str]],
    logic_types: Mapping[str, str],
    ignore_rv0678_neutral_fs: bool = False,
) -> list[dict[str, str]]:
    augmented: list[dict[str, str]] = []
    for row in rows:
        output = dict(row)
        output["original_current_evaluator_scope"] = row.get(
            "current_evaluator_scope", ""
        )
        output["original_amr_hunter_prediction"] = row.get("amr_hunter_prediction", "")
        output["original_matched_status"] = row.get("matched_status", "")
        output["original_correct"] = row.get("correct", "")
        output["original_gap_reason"] = row.get("gap_reason", "")
        output["call_source"] = (
            "direct_record_match"
            if row.get("matched_status") == "matched_final"
            else ""
        )
        output["consequence_inference"] = ""
        output["inferred_functional_state"] = ""
        output["inferred_logic_type"] = ""
        if can_infer_lof(row, logic_types):
            logic_type = logic_types[str(row.get("gene") or "")]
            decision = resolve_consequence_level_decision(
                {**row, "logic_type": logic_type},
                drug="Bedaquiline",
                ignore_rv0678_neutral_fs=ignore_rv0678_neutral_fs,
            )
            if decision and decision.resistance_phenotype in {"RESISTANT", "SENSITIVE"}:
                prediction = (
                    "R" if decision.resistance_phenotype == "RESISTANT" else "S"
                )
                original_evidence = str(row.get("evidence_code") or "")
                evidence_codes = [
                    code for code in (original_evidence, decision.evidence_code) if code
                ]
                output["current_evaluator_scope"] = "true"
                output["amr_hunter_prediction"] = prediction
                output["matched_status"] = "consequence_inferred"
                output["correct"] = (
                    "true" if prediction == row.get("who_truth") else "false"
                )
                output["gap_reason"] = ""
                output["record_count"] = row.get("record_count", "0") or "0"
                output["source_statuses"] = row.get("source_statuses", "")
                output["evidence_code"] = "|".join(dict.fromkeys(evidence_codes))
                output["source_records"] = (
                    f"no_direct_final_rs_record;inferred_from_{row.get('mutation_class')}_as_loss_of_function;logic_type={logic_type}"
                )
                output["call_source"] = "catalogue_consequence_inferred"
                output["consequence_inference"] = "loss_of_function"
                output["inferred_functional_state"] = decision.functional_state
                output["inferred_logic_type"] = logic_type
        augmented.append(output)
    return augmented


def is_primary_evaluable(row: Mapping[str, str]) -> bool:
    return row.get("who_truth") in {"R", "S"} and row.get("benchmark_role") == "primary"


def is_current_scope(row: Mapping[str, str]) -> bool:
    return is_primary_evaluable(row) and row.get("current_evaluator_scope") == "true"


def is_matched(row: Mapping[str, str]) -> bool:
    return is_primary_evaluable(row) and row.get("matched_status") in MATCHED_STATUSES


def compact_counter(counter: Counter[str], limit: int = 8) -> str:
    return ";".join((f"{key}:{value}" for key, value in counter.most_common(limit)))


def summarize_group(
    scope_type: str, scope_value: str, rows: Sequence[Mapping[str, str]]
) -> dict[str, object]:
    primary = [row for row in rows if is_primary_evaluable(row)]
    current_scope = [row for row in rows if is_current_scope(row)]
    matched = [row for row in rows if is_matched(row)]
    direct = [row for row in matched if row.get("call_source") == "direct_record_match"]
    inferred = [
        row
        for row in matched
        if row.get("call_source") == "catalogue_consequence_inferred"
    ]
    correct = [row for row in matched if row.get("correct") == "true"]
    wrong = [row for row in matched if row.get("correct") == "false"]
    unmatched = [row for row in primary if not is_matched(row)]
    cm = Counter(
        ((row.get("who_truth"), row.get("amr_hunter_prediction")) for row in matched)
    )
    gap_counts = Counter(
        (row.get("gap_reason") for row in unmatched if row.get("gap_reason"))
    )
    return {
        "scope_type": scope_type,
        "scope_value": scope_value,
        "who_rows": len(rows),
        "primary_evaluable_rows": len(primary),
        "current_evaluator_scope_rows": len(current_scope),
        "matched_rows": len(matched),
        "direct_matched_rows": len(direct),
        "lof_inferred_rows": len(inferred),
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
    rows: Sequence[Mapping[str, str]], keys: Sequence[str]
) -> dict[tuple[str, ...], list[Mapping[str, str]]]:
    groups: dict[tuple[str, ...], list[Mapping[str, str]]] = defaultdict(list)
    for row in rows:
        groups[tuple((row.get(key) or "BLANK" for key in keys))].append(row)
    return groups


def build_stratified_rows(rows: Sequence[Mapping[str, str]]) -> list[dict[str, object]]:
    output: list[dict[str, object]] = [summarize_group("overall_primary", "all", rows)]
    for key, items in sorted(group_by(rows, ["gene"]).items()):
        output.append(summarize_group("gene", key[0], items))
    primary_rows = [row for row in rows if row.get("benchmark_role") == "primary"]
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
        row for row in rows if is_primary_evaluable(row) and (not is_matched(row))
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
                    (row.get("who_row_id", "") for row in examples)
                ),
                "example_mutations": "|".join(
                    (row.get("mutation", "") for row in examples)
                ),
            }
        )
    return output


def main() -> int:
    args = parse_args()
    rows, input_fields = read_rows(args.input)
    logic_types = load_logic_types(args.config)
    augmented = augment_rows(
        rows, logic_types, ignore_rv0678_neutral_fs=args.ignore_rv0678_neutral_fs
    )
    summary = [summarize_group("overall_primary", "all", augmented)]
    stratified = build_stratified_rows(augmented)
    gaps = build_gap_rows(augmented)
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
        args.output_dir / "who_bdq_per_row_benchmark_lof_augmented.csv",
        augmented,
        per_row_fields,
    )
    write_csv(
        args.output_dir / "who_bdq_benchmark_summary_lof_augmented.csv",
        summary,
        summary_fields,
    )
    write_csv(
        args.output_dir / "who_bdq_stratified_benchmark_lof_augmented.csv",
        stratified,
        summary_fields,
    )
    write_csv(
        args.output_dir / "who_bdq_gap_breakdown_lof_augmented.csv", gaps, gap_fields
    )
    overall = summary[0]
    print(f"per_row={args.output_dir / 'who_bdq_per_row_benchmark_lof_augmented.csv'}")
    print(f"summary={args.output_dir / 'who_bdq_benchmark_summary_lof_augmented.csv'}")
    print(
        "matched={matched_rows} direct={direct_matched_rows} lof_inferred={lof_inferred_rows} coverage={coverage} accuracy={matched_call_accuracy}".format(
            **overall
        )
    )
    print(f"remaining_unmatched={overall['unmatched_rows']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
