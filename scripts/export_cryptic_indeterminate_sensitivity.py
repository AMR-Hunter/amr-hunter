from __future__ import annotations
import argparse
import csv
from collections import defaultdict
from pathlib import Path

DEFAULT_ISOLATES = Path(
    "reports/publication/external_bdq_phenotype_validation_isolates.csv"
)
DEFAULT_VARIANTS = Path(
    "reports/publication/external_bdq_phenotype_validation_variant_calls.csv"
)
DEFAULT_SUMMARY_OUTPUT = Path(
    "reports/publication/cryptic_indeterminate_sensitivity_summary.csv"
)
DEFAULT_ISOLATE_OUTPUT = Path(
    "reports/publication/cryptic_indeterminate_sensitivity_isolates.csv"
)
CORE_GENES = {"Rv0678", "atpE", "pepQ"}
SECONDARY_GENES = {"glpK", "Rv1979c", "mmpL5", "mmpS5"}
CONTEXT_GENES = {"mtrA", "mtrB"}
MTRB_CONTEXT_VARIANTS = {("mtrB", "p.Met517Leu"), ("mtrB", "p.Pro18Ser")}
SCENARIOS = [
    {
        "scenario": "baseline",
        "description": "No indeterminate override; use calibrated CRyPTIC validation calls as exported.",
        "genes": set(),
        "variant_keys": set(),
    },
    {
        "scenario": "core_unsupported",
        "description": "Current S calls become INDETERMINATE when an unsupported core BDQ gene variant is present.",
        "genes": CORE_GENES,
        "variant_keys": set(),
    },
    {
        "scenario": "core_secondary_unsupported",
        "description": "Current S calls become INDETERMINATE when unsupported core or secondary BDQ gene variants are present.",
        "genes": CORE_GENES | SECONDARY_GENES,
        "variant_keys": set(),
    },
    {
        "scenario": "mtrB_context",
        "description": "Current S calls become INDETERMINATE when unsupported mtrB Met517Leu or Pro18Ser is present.",
        "genes": set(),
        "variant_keys": MTRB_CONTEXT_VARIANTS,
    },
    {
        "scenario": "core_secondary_mtrB_context",
        "description": "Current S calls become INDETERMINATE for unsupported core/secondary variants or mtrB context variants.",
        "genes": CORE_GENES | SECONDARY_GENES,
        "variant_keys": MTRB_CONTEXT_VARIANTS,
    },
    {
        "scenario": "any_unsupported_bdq_scope",
        "description": "Current S calls become INDETERMINATE when any unsupported configured BDQ-scope variant is present.",
        "genes": CORE_GENES | SECONDARY_GENES | CONTEXT_GENES,
        "variant_keys": set(),
    },
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export indeterminate sensitivity scenarios for CRyPTIC BDQ validation."
    )
    parser.add_argument("--isolates", type=Path, default=DEFAULT_ISOLATES)
    parser.add_argument("--variants", type=Path, default=DEFAULT_VARIANTS)
    parser.add_argument("--summary-output", type=Path, default=DEFAULT_SUMMARY_OUTPUT)
    parser.add_argument("--isolate-output", type=Path, default=DEFAULT_ISOLATE_OUTPUT)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def has_unsupported_trigger(
    variants: list[dict[str, str]], genes: set[str], variant_keys: set[tuple[str, str]]
) -> tuple[bool, list[str]]:
    labels: list[str] = []
    for row in variants:
        if row["call_source"] != "unsupported_variant":
            continue
        key = (row["gene"], row["mutation"])
        if row["gene"] in genes or key in variant_keys:
            labels.append(f"{row['gene']}:{row['mutation']}")
    return (bool(labels), sorted(labels))


def scenario_prediction(
    isolate: dict[str, str], variants: list[dict[str, str]], scenario: dict[str, object]
) -> tuple[str, str]:
    baseline = isolate["amr_hunter_isolate_prediction"]
    if scenario["scenario"] == "baseline" or baseline != "S":
        return (baseline, "")
    triggered, labels = has_unsupported_trigger(
        variants, scenario["genes"], scenario["variant_keys"]
    )
    if not triggered:
        return (baseline, "")
    return ("INDETERMINATE", ";".join(labels))


def summarize_predictions(
    isolates: list[dict[str, str]],
    variants_by_isolate: dict[str, list[dict[str, str]]],
    scenario: dict[str, object],
) -> tuple[dict[str, str], list[dict[str, str]]]:
    total = len(isolates)
    phenotype_r = sum((1 for row in isolates if row["bdq_binary_phenotype"] == "R"))
    phenotype_s = sum((1 for row in isolates if row["bdq_binary_phenotype"] == "S"))
    tp = tn = fp = fn = 0
    pred_r = pred_s = pred_unknown = pred_indeterminate = 0
    r_indeterminate = s_indeterminate = 0
    isolate_rows: list[dict[str, str]] = []
    for isolate in isolates:
        prediction, trigger_labels = scenario_prediction(
            isolate, variants_by_isolate.get(isolate["unique_id"], []), scenario
        )
        phenotype = isolate["bdq_binary_phenotype"]
        if prediction == "R":
            pred_r += 1
            if phenotype == "R":
                tp += 1
            elif phenotype == "S":
                fp += 1
        elif prediction == "S":
            pred_s += 1
            if phenotype == "R":
                fn += 1
            elif phenotype == "S":
                tn += 1
        elif prediction == "INDETERMINATE":
            pred_indeterminate += 1
            if phenotype == "R":
                r_indeterminate += 1
            elif phenotype == "S":
                s_indeterminate += 1
        else:
            pred_unknown += 1
        if prediction != isolate["amr_hunter_isolate_prediction"]:
            isolate_rows.append(
                {
                    "scenario": str(scenario["scenario"]),
                    "unique_id": isolate["unique_id"],
                    "ena_sample": isolate["ena_sample"],
                    "bdq_binary_phenotype": phenotype,
                    "bdq_mic": isolate["bdq_mic"],
                    "baseline_prediction": isolate["amr_hunter_isolate_prediction"],
                    "scenario_prediction": prediction,
                    "trigger_variants": trigger_labels,
                    "baseline_call_source_summary": isolate["call_source_summary"],
                }
            )
    hard_calls = tp + tn + fp + fn
    unknown_like = pred_unknown + pred_indeterminate

    def ratio(num: int, den: int) -> str:
        return f"{num / den:.4f}" if den else ""

    summary = {
        "scenario": str(scenario["scenario"]),
        "description": str(scenario["description"]),
        "total_isolates": str(total),
        "phenotype_r": str(phenotype_r),
        "phenotype_s": str(phenotype_s),
        "hard_call_isolates": str(hard_calls),
        "unknown_or_indeterminate_isolates": str(unknown_like),
        "coverage": ratio(hard_calls, total),
        "hard_call_accuracy": ratio(tp + tn, hard_calls),
        "tp": str(tp),
        "tn": str(tn),
        "fp": str(fp),
        "fn": str(fn),
        "pred_R": str(pred_r),
        "pred_S": str(pred_s),
        "pred_UNKNOWN": str(pred_unknown),
        "pred_INDETERMINATE": str(pred_indeterminate),
        "r_indeterminate": str(r_indeterminate),
        "s_indeterminate": str(s_indeterminate),
        "hard_call_sensitivity": ratio(tp, tp + fn),
        "hard_call_specificity": ratio(tn, tn + fp),
        "r_not_called_sensitive": ratio(tp + r_indeterminate, phenotype_r),
        "s_not_called_resistant": ratio(tn + s_indeterminate, phenotype_s),
        "fn_reduction_vs_baseline": "",
        "fp_reduction_vs_baseline": "",
        "coverage_delta_vs_baseline": "",
    }
    return (summary, isolate_rows)


def add_baseline_deltas(rows: list[dict[str, str]]) -> None:
    baseline = next((row for row in rows if row["scenario"] == "baseline"), None)
    if baseline is None:
        return
    baseline_fn = int(baseline["fn"])
    baseline_fp = int(baseline["fp"])
    baseline_coverage = float(baseline["coverage"])
    for row in rows:
        row["fn_reduction_vs_baseline"] = str(baseline_fn - int(row["fn"]))
        row["fp_reduction_vs_baseline"] = str(baseline_fp - int(row["fp"]))
        coverage = float(row["coverage"]) if row["coverage"] else 0.0
        row["coverage_delta_vs_baseline"] = f"{coverage - baseline_coverage:.4f}"


def main() -> int:
    args = parse_args()
    isolates = read_csv(args.isolates)
    variants = read_csv(args.variants)
    variants_by_isolate: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in variants:
        variants_by_isolate[row["unique_id"]].append(row)
    summary_rows: list[dict[str, str]] = []
    isolate_rows: list[dict[str, str]] = []
    for scenario in SCENARIOS:
        summary, changed_isolates = summarize_predictions(
            isolates, variants_by_isolate, scenario
        )
        summary_rows.append(summary)
        isolate_rows.extend(changed_isolates)
    add_baseline_deltas(summary_rows)
    write_csv(args.summary_output, summary_rows)
    write_csv(args.isolate_output, isolate_rows)
    print(f"summary={args.summary_output}")
    print(f"summary_rows={len(summary_rows)}")
    print(f"isolates={args.isolate_output}")
    print(f"changed_isolate_rows={len(isolate_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
