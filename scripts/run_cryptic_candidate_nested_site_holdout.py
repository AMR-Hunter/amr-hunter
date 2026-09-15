from __future__ import annotations
import argparse
import re
import sys
from collections import Counter
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
from scripts.export_cryptic_candidate_nested_validation import (
    NESTED_MODEL,
    REFERENCE_MODEL,
    add_candidate_columns,
    select_train_candidates,
    unsupported_variants_by_isolate,
)
from scripts.export_cryptic_supervised_upper_bound import (
    TARGET_SPECIFICITIES,
    build_feature_table,
    choose_threshold,
    confusion,
    read_csv,
    ratio,
    score_row,
    selected_features,
    train_bernoulli_nb,
    write_csv,
)

DEFAULT_ISOLATES = Path(
    "reports/publication/external_bdq_phenotype_validation_isolates.csv"
)
DEFAULT_VARIANTS = Path(
    "reports/publication/external_bdq_phenotype_validation_variant_calls.csv"
)
DEFAULT_SUMMARY_OUTPUT = Path(
    "reports/publication/cryptic_candidate_nested_site_holdout_summary.csv"
)
DEFAULT_METRICS_OUTPUT = Path(
    "reports/publication/cryptic_candidate_nested_site_holdout_metrics.csv"
)
DEFAULT_PREDICTION_OUTPUT = Path(
    "reports/publication/cryptic_candidate_nested_site_holdout_predictions.csv"
)
DEFAULT_SELECTION_OUTPUT = Path(
    "reports/publication/cryptic_candidate_nested_site_holdout_selection.csv"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Evaluate candidate unsupported variants with leave-one-site-out nested selection."
    )
    parser.add_argument("--isolates", type=Path, default=DEFAULT_ISOLATES)
    parser.add_argument("--variants", type=Path, default=DEFAULT_VARIANTS)
    parser.add_argument("--summary-output", type=Path, default=DEFAULT_SUMMARY_OUTPUT)
    parser.add_argument("--metrics-output", type=Path, default=DEFAULT_METRICS_OUTPUT)
    parser.add_argument(
        "--prediction-output", type=Path, default=DEFAULT_PREDICTION_OUTPUT
    )
    parser.add_argument(
        "--selection-output", type=Path, default=DEFAULT_SELECTION_OUTPUT
    )
    parser.add_argument("--min-train-r-carriers", type=int, default=2)
    parser.add_argument("--min-train-fn-carriers", type=int, default=1)
    parser.add_argument("--min-r-fraction", type=float, default=0.05)
    parser.add_argument("--min-odds-ratio", type=float, default=2.0)
    parser.add_argument("--max-p-value", type=float, default=0.05)
    parser.add_argument("--min-train-positive", type=int, default=1)
    return parser.parse_args()


def infer_site(unique_id: str) -> str:
    match = re.match("(site\\.\\d+)", unique_id)
    return match.group(1) if match else "unknown"


def add_site(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    out = []
    for row in rows:
        copied = dict(row)
        copied["site"] = infer_site(row["unique_id"])
        out.append(copied)
    return out


def group_counts(rows: list[dict[str, str]], indices: list[int]) -> Counter[str]:
    counts: Counter[str] = Counter()
    for idx in indices:
        counts[rows[idx]["bdq_binary_phenotype"]] += 1
    return counts


def summarize_counts(
    model: str, target: float, counts: Counter[str], groups: int
) -> dict[str, str]:
    total = counts["tp"] + counts["tn"] + counts["fp"] + counts["fn"]
    return {
        "model": model,
        "target_specificity": f"{target:.2f}",
        "heldout_groups": str(groups),
        "total": str(total),
        "tp": str(counts["tp"]),
        "tn": str(counts["tn"]),
        "fp": str(counts["fp"]),
        "fn": str(counts["fn"]),
        "recall": ratio(counts["tp"], counts["tp"] + counts["fn"]),
        "specificity": ratio(counts["tn"], counts["tn"] + counts["fp"]),
        "precision": ratio(counts["tp"], counts["tp"] + counts["fp"]),
        "accuracy": ratio(counts["tp"] + counts["tn"], total),
    }


def run_site_holdout_nested(
    rows: list[dict[str, str]],
    base_feature_names: list[str],
    unsupported_by_isolate: dict[str, set[str]],
    min_train_r_carriers: int,
    min_train_fn_carriers: int,
    min_r_fraction: float,
    min_odds_ratio: float,
    max_p_value: float,
    min_train_positive: int,
) -> tuple[
    list[dict[str, str]],
    list[dict[str, str]],
    list[dict[str, str]],
    list[dict[str, str]],
]:
    sites = sorted({row["site"] for row in rows})
    site_to_indices = {
        site: [idx for idx, row in enumerate(rows) if row["site"] == site]
        for site in sites
    }
    reference_features = selected_features(
        base_feature_names, exclude_groups={"secondary"}
    )
    aggregate = {
        REFERENCE_MODEL: {target: Counter() for target in TARGET_SPECIFICITIES},
        NESTED_MODEL: {target: Counter() for target in TARGET_SPECIFICITIES},
    }
    evaluated_groups = {REFERENCE_MODEL: Counter(), NESTED_MODEL: Counter()}
    metric_rows: list[dict[str, str]] = []
    prediction_rows: list[dict[str, str]] = []
    selection_rows: list[dict[str, str]] = []
    for site in sites:
        test_indices = site_to_indices[site]
        test_set = set(test_indices)
        train_indices = [idx for idx in range(len(rows)) if idx not in test_set]
        train_counts = group_counts(rows, train_indices)
        test_counts = group_counts(rows, test_indices)
        if train_counts["R"] < min_train_positive:
            continue
        selected_candidates = select_train_candidates(
            rows,
            train_indices,
            unsupported_by_isolate,
            min_train_r_carriers,
            min_train_fn_carriers,
            min_r_fraction,
            min_odds_ratio,
            max_p_value,
        )
        candidate_variants = [row["candidate"] for row in selected_candidates]
        candidate_features = add_candidate_columns(
            rows, candidate_variants, unsupported_by_isolate
        )
        for rank, row in enumerate(selected_candidates, start=1):
            selection_rows.append({"heldout_site": site, "rank": str(rank), **row})
        model_features = {
            REFERENCE_MODEL: reference_features,
            NESTED_MODEL: [*reference_features, *candidate_features],
        }
        for model_name, features in model_features.items():
            model = train_bernoulli_nb(rows, train_indices, features)
            scores = {
                idx: score_row(rows[idx], model, features)
                for idx in train_indices + test_indices
            }
            thresholds = {
                target: choose_threshold(rows, train_indices, scores, target)
                for target in TARGET_SPECIFICITIES
            }
            for target, threshold in thresholds.items():
                counts = confusion(rows, test_indices, scores, threshold)
                aggregate[model_name][target].update(counts)
                evaluated_groups[model_name][target] += 1
                metric_row = summarize_counts(model_name, target, counts, groups=1)
                metric_row.update(
                    {
                        "heldout_site": site,
                        "site_total": str(len(test_indices)),
                        "site_r": str(test_counts["R"]),
                        "site_s": str(test_counts["S"]),
                        "train_r": str(train_counts["R"]),
                        "train_s": str(train_counts["S"]),
                        "selected_candidate_count": str(len(candidate_features)),
                        "selected_candidates": ";".join(candidate_variants),
                        "threshold": f"{threshold:.6f}",
                    }
                )
                metric_rows.append(metric_row)
            for idx in test_indices:
                pred_row = {
                    "model": model_name,
                    "heldout_site": site,
                    "unique_id": rows[idx]["unique_id"],
                    "bdq_binary_phenotype": rows[idx]["bdq_binary_phenotype"],
                    "score": f"{scores[idx]:.6f}",
                    "selected_candidate_count": str(len(candidate_features)),
                    "selected_candidates": ";".join(candidate_variants),
                }
                for target, threshold in thresholds.items():
                    pred_row[f"pred_at_spec_{target:.2f}"] = (
                        "R" if scores[idx] >= threshold else "S"
                    )
                prediction_rows.append(pred_row)
    summary_rows = []
    for model_name in (REFERENCE_MODEL, NESTED_MODEL):
        for target in TARGET_SPECIFICITIES:
            summary_rows.append(
                summarize_counts(
                    model_name,
                    target,
                    aggregate[model_name][target],
                    groups=evaluated_groups[model_name][target],
                )
            )
    return (summary_rows, metric_rows, prediction_rows, selection_rows)


def main() -> int:
    args = parse_args()
    isolates = read_csv(args.isolates)
    variants = read_csv(args.variants)
    rows, base_feature_names = build_feature_table(isolates, variants)
    rows = add_site(rows)
    unsupported = unsupported_variants_by_isolate(variants)
    summary_rows, metric_rows, prediction_rows, selection_rows = (
        run_site_holdout_nested(
            rows,
            base_feature_names,
            unsupported,
            args.min_train_r_carriers,
            args.min_train_fn_carriers,
            args.min_r_fraction,
            args.min_odds_ratio,
            args.max_p_value,
            args.min_train_positive,
        )
    )
    write_csv(args.summary_output, summary_rows)
    write_csv(args.metrics_output, metric_rows)
    write_csv(args.prediction_output, prediction_rows)
    write_csv(args.selection_output, selection_rows)
    print(f"summary={args.summary_output}")
    print(f"metrics={args.metrics_output}")
    print(f"predictions={args.prediction_output}")
    print(f"prediction_rows={len(prediction_rows)}")
    print(f"selection={args.selection_output}")
    print(f"selection_rows={len(selection_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
