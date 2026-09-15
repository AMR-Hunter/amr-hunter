from __future__ import annotations
import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
from scripts.export_cryptic_supervised_upper_bound import (
    MODEL_SPECS,
    TARGET_SPECIFICITIES,
    choose_threshold,
    confusion,
    recall,
    score_row,
    selected_features,
    specificity,
    train_bernoulli_nb,
)

DEFAULT_FEATURES = Path("reports/publication/cryptic_isolate_feature_table.csv")
DEFAULT_SUMMARY_OUTPUT = Path(
    "reports/publication/cryptic_supervised_site_holdout_summary.csv"
)
DEFAULT_METRICS_OUTPUT = Path(
    "reports/publication/cryptic_supervised_site_holdout_metrics.csv"
)
DEFAULT_PREDICTION_OUTPUT = Path(
    "reports/publication/cryptic_supervised_site_holdout_predictions.csv"
)
METADATA_COLUMNS = {
    "unique_id",
    "ena_sample",
    "bdq_binary_phenotype",
    "label",
    "bdq_mic",
    "baseline_prediction",
    "site",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Evaluate supervised CRyPTIC BDQ models with leave-one-site-out validation."
    )
    parser.add_argument("--features", type=Path, default=DEFAULT_FEATURES)
    parser.add_argument("--summary-output", type=Path, default=DEFAULT_SUMMARY_OUTPUT)
    parser.add_argument("--metrics-output", type=Path, default=DEFAULT_METRICS_OUTPUT)
    parser.add_argument(
        "--prediction-output", type=Path, default=DEFAULT_PREDICTION_OUTPUT
    )
    parser.add_argument("--min-train-positive", type=int, default=1)
    parser.add_argument("--min-test-total", type=int, default=1)
    return parser.parse_args()


def read_csv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return (list(reader), list(reader.fieldnames or []))


def write_csv(
    path: Path, rows: list[dict[str, str]], fieldnames: list[str] | None = None
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


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


def feature_names(fieldnames: list[str]) -> list[str]:
    return [name for name in fieldnames if name not in METADATA_COLUMNS]


def ratio(num: int, den: int) -> str:
    return f"{num / den:.4f}" if den else ""


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


def group_counts(rows: list[dict[str, str]], indices: list[int]) -> Counter[str]:
    counts: Counter[str] = Counter()
    for idx in indices:
        counts[rows[idx]["bdq_binary_phenotype"]] += 1
    return counts


def run_site_holdout(
    rows: list[dict[str, str]],
    all_features: list[str],
    min_train_positive: int = 1,
    min_test_total: int = 1,
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    sites = sorted({row["site"] for row in rows})
    site_to_indices = {
        site: [idx for idx, row in enumerate(rows) if row["site"] == site]
        for site in sites
    }
    summary_rows: list[dict[str, str]] = []
    metric_rows: list[dict[str, str]] = []
    prediction_rows: list[dict[str, str]] = []
    for model_name, spec in MODEL_SPECS.items():
        features = selected_features(
            all_features, include_groups=spec["include"], exclude_groups=spec["exclude"]
        )
        aggregate = {target: Counter() for target in TARGET_SPECIFICITIES}
        evaluated_groups = Counter()
        for site in sites:
            test_indices = site_to_indices[site]
            if len(test_indices) < min_test_total:
                continue
            train_indices = [
                idx for idx in range(len(rows)) if idx not in set(test_indices)
            ]
            train_counts = group_counts(rows, train_indices)
            test_counts = group_counts(rows, test_indices)
            if train_counts["R"] < min_train_positive:
                continue
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
                aggregate[target].update(counts)
                evaluated_groups[target] += 1
                metric_row = summarize_counts(model_name, target, counts, groups=1)
                metric_row.update(
                    {
                        "heldout_site": site,
                        "site_total": str(len(test_indices)),
                        "site_r": str(test_counts["R"]),
                        "site_s": str(test_counts["S"]),
                        "train_r": str(train_counts["R"]),
                        "train_s": str(train_counts["S"]),
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
                }
                for target, threshold in thresholds.items():
                    pred_row[f"pred_at_spec_{target:.2f}"] = (
                        "R" if scores[idx] >= threshold else "S"
                    )
                prediction_rows.append(pred_row)
        for target in TARGET_SPECIFICITIES:
            summary_rows.append(
                summarize_counts(
                    model_name,
                    target,
                    aggregate[target],
                    groups=evaluated_groups[target],
                )
            )
    return (summary_rows, metric_rows, prediction_rows)


def main() -> int:
    args = parse_args()
    rows, fieldnames = read_csv(args.features)
    rows = add_site(rows)
    features = feature_names(fieldnames)
    summary_rows, metric_rows, prediction_rows = run_site_holdout(
        rows,
        features,
        min_train_positive=args.min_train_positive,
        min_test_total=args.min_test_total,
    )
    write_csv(args.summary_output, summary_rows)
    write_csv(args.metrics_output, metric_rows)
    write_csv(args.prediction_output, prediction_rows)
    print(f"summary={args.summary_output}")
    print(f"metrics={args.metrics_output}")
    print(f"predictions={args.prediction_output}")
    print(f"prediction_rows={len(prediction_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
