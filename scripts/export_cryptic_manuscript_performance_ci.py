from __future__ import annotations
import argparse
import csv
import random
from collections import defaultdict
from pathlib import Path

DEFAULT_TABLE = Path("reports/publication/cryptic_manuscript_performance_table.csv")
DEFAULT_ISOLATES = Path(
    "reports/publication/external_bdq_phenotype_validation_isolates.csv"
)
DEFAULT_FIVE_LEVEL = Path("reports/publication/cryptic_five_level_isolates.csv")
DEFAULT_CV_PREDICTIONS = Path(
    "reports/publication/cryptic_supervised_cv_predictions.csv"
)
DEFAULT_SITE_PREDICTIONS = Path(
    "reports/publication/cryptic_supervised_site_holdout_predictions.csv"
)
DEFAULT_CI_OUTPUT = Path(
    "reports/publication/cryptic_manuscript_performance_bootstrap_ci.csv"
)
DEFAULT_TABLE_CI_OUTPUT = Path(
    "reports/publication/cryptic_manuscript_performance_table_with_ci.csv"
)
DEFAULT_TABLE_CI_MD_OUTPUT = Path(
    "reports/publication/cryptic_manuscript_performance_table_with_ci.md"
)
METRICS = ("recall", "specificity", "precision")
CI_FIELDNAMES = [
    "analysis",
    "model",
    "threshold_target",
    "metric",
    "estimate",
    "bootstrap_ci_low",
    "bootstrap_ci_high",
    "bootstrap_iterations",
    "confidence_level",
    "resampling_unit",
    "resampled_units",
    "notes",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export bootstrap CIs for CRyPTIC manuscript performance metrics."
    )
    parser.add_argument("--table", type=Path, default=DEFAULT_TABLE)
    parser.add_argument("--isolates", type=Path, default=DEFAULT_ISOLATES)
    parser.add_argument("--five-level", type=Path, default=DEFAULT_FIVE_LEVEL)
    parser.add_argument("--cv-predictions", type=Path, default=DEFAULT_CV_PREDICTIONS)
    parser.add_argument(
        "--site-predictions", type=Path, default=DEFAULT_SITE_PREDICTIONS
    )
    parser.add_argument("--ci-output", type=Path, default=DEFAULT_CI_OUTPUT)
    parser.add_argument("--table-ci-output", type=Path, default=DEFAULT_TABLE_CI_OUTPUT)
    parser.add_argument(
        "--table-ci-md-output", type=Path, default=DEFAULT_TABLE_CI_MD_OUTPUT
    )
    parser.add_argument("--iterations", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=271828)
    parser.add_argument("--confidence", type=float, default=0.95)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def is_positive(value: str) -> bool:
    return value == "R"


def counts_for(truth: str, pred: str | None) -> tuple[int, int, int, int]:
    if pred is None:
        return (0, 0, 0, 0)
    truth_r = is_positive(truth)
    pred_r = is_positive(pred)
    if truth_r and pred_r:
        return (1, 0, 0, 0)
    if not truth_r and (not pred_r):
        return (0, 1, 0, 0)
    if not truth_r and pred_r:
        return (0, 0, 1, 0)
    return (0, 0, 0, 1)


def add_counts(
    a: tuple[int, int, int, int], b: tuple[int, int, int, int]
) -> tuple[int, int, int, int]:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3])


def metric_values(counts: tuple[int, int, int, int]) -> dict[str, float | None]:
    tp, tn, fp, fn = counts
    return {
        "recall": tp / (tp + fn) if tp + fn else None,
        "specificity": tn / (tn + fp) if tn + fp else None,
        "precision": tp / (tp + fp) if tp + fp else None,
    }


def format_float(value: float | None) -> str:
    return f"{value:.4f}" if value is not None else ""


def percentile(values: list[float], quantile: float) -> float:
    if not values:
        raise ValueError("Cannot compute percentile of an empty list")
    ordered = sorted(values)
    position = quantile * (len(ordered) - 1)
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def bootstrap_metric_cis(
    grouped_counts: list[tuple[int, int, int, int]],
    iterations: int,
    confidence: float,
    rng: random.Random,
) -> dict[str, tuple[float | None, float | None]]:
    alpha = 1 - confidence
    n = len(grouped_counts)
    sampled_metrics: dict[str, list[float]] = {metric: [] for metric in METRICS}
    for _ in range(iterations):
        total = (0, 0, 0, 0)
        for _ in range(n):
            total = add_counts(total, grouped_counts[rng.randrange(n)])
        values = metric_values(total)
        for metric, value in values.items():
            if value is not None:
                sampled_metrics[metric].append(value)
    return {
        metric: (
            percentile(values, alpha / 2) if values else None,
            percentile(values, 1 - alpha / 2) if values else None,
        )
        for metric, values in sampled_metrics.items()
    }


def frozen_groups(rows: list[dict[str, str]]) -> list[tuple[int, int, int, int]]:
    groups = []
    for row in rows:
        prediction = row["amr_hunter_isolate_prediction"]
        pred = prediction if prediction in {"R", "S"} else None
        groups.append(counts_for(row["bdq_binary_phenotype"], pred))
    return groups


def five_level_groups(
    rows: list[dict[str, str]], positive_labels: set[str]
) -> list[tuple[int, int, int, int]]:
    groups = []
    for row in rows:
        pred = "R" if row["five_level_label"] in positive_labels else "S"
        groups.append(counts_for(row["bdq_binary_phenotype"], pred))
    return groups


def prediction_groups(
    rows: list[dict[str, str]], model: str, target: str, cluster_by_isolate: bool
) -> list[tuple[int, int, int, int]]:
    pred_column = f"pred_at_spec_{target}"
    selected = [row for row in rows if row["model"] == model]
    if cluster_by_isolate:
        grouped: dict[str, tuple[int, int, int, int]] = defaultdict(
            lambda: (0, 0, 0, 0)
        )
        for row in selected:
            grouped[row["unique_id"]] = add_counts(
                grouped[row["unique_id"]],
                counts_for(row["bdq_binary_phenotype"], row[pred_column]),
            )
        return list(grouped.values())
    return [
        counts_for(row["bdq_binary_phenotype"], row[pred_column]) for row in selected
    ]


def groups_for_analysis(
    row: dict[str, str],
    isolates: list[dict[str, str]],
    five_level: list[dict[str, str]],
    cv_predictions: list[dict[str, str]],
    site_predictions: list[dict[str, str]],
) -> tuple[list[tuple[int, int, int, int]], str, str]:
    analysis = row["analysis"]
    if analysis == "Frozen binary AMR-Hunter":
        return (
            frozen_groups(isolates),
            "isolate",
            "UNKNOWN/no-call isolates contribute no TP/TN/FP/FN; coverage is reported separately.",
        )
    if analysis == "Five-level hard R only":
        return (
            five_level_groups(five_level, {"R"}),
            "isolate",
            "Indeterminate and Likely R labels are treated as non-resistant in this binary projection.",
        )
    if analysis == "Five-level R or Likely R":
        return (
            five_level_groups(five_level, {"R", "Likely R"}),
            "isolate",
            "Indeterminate labels are treated as non-resistant in this binary projection.",
        )
    if analysis.startswith("Repeated 5x5 CV"):
        return (
            prediction_groups(
                cv_predictions,
                row["model"],
                row["threshold_target"],
                cluster_by_isolate=True,
            ),
            "isolate_cluster_with_5_cv_predictions",
            "Cluster bootstrap resamples isolates and keeps their repeated CV predictions together.",
        )
    if analysis.startswith("Leave-one-site-out"):
        return (
            prediction_groups(
                site_predictions,
                row["model"],
                row["threshold_target"],
                cluster_by_isolate=True,
            ),
            "isolate",
            "Bootstrap resamples held-out isolate predictions.",
        )
    raise ValueError(f"Unsupported analysis row: {analysis}")


def build_ci_rows(
    table_rows: list[dict[str, str]],
    isolates: list[dict[str, str]],
    five_level: list[dict[str, str]],
    cv_predictions: list[dict[str, str]],
    site_predictions: list[dict[str, str]],
    iterations: int,
    confidence: float,
    seed: int,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    rng = random.Random(seed)
    ci_rows = []
    table_with_ci = []
    for row in table_rows:
        groups, unit, notes = groups_for_analysis(
            row, isolates, five_level, cv_predictions, site_predictions
        )
        estimates = metric_values(
            tuple((int(row[name]) for name in ("tp", "tn", "fp", "fn")))
        )
        cis = bootstrap_metric_cis(groups, iterations, confidence, rng)
        enriched = dict(row)
        for metric in METRICS:
            low, high = cis[metric]
            ci_rows.append(
                {
                    "analysis": row["analysis"],
                    "model": row["model"],
                    "threshold_target": row["threshold_target"],
                    "metric": metric,
                    "estimate": format_float(estimates[metric]),
                    "bootstrap_ci_low": format_float(low),
                    "bootstrap_ci_high": format_float(high),
                    "bootstrap_iterations": str(iterations),
                    "confidence_level": f"{confidence:.2f}",
                    "resampling_unit": unit,
                    "resampled_units": str(len(groups)),
                    "notes": notes,
                }
            )
            enriched[f"{metric}_ci"] = (
                f"{format_float(estimates[metric])} ({format_float(low)}-{format_float(high)})"
                if estimates[metric] is not None
                and low is not None
                and (high is not None)
                else ""
            )
        table_with_ci.append(enriched)
    return (ci_rows, table_with_ci)


def markdown_table(rows: list[dict[str, str]]) -> str:
    columns = [
        "analysis",
        "evaluation_design",
        "model",
        "threshold_target",
        "recall_ci",
        "specificity_ci",
        "precision_ci",
        "coverage",
        "indeterminate",
    ]
    lines = [
        "# CRyPTIC BDQ Manuscript Performance Table With Bootstrap CIs",
        "",
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for row in rows:
        lines.append(
            "| " + " | ".join((row.get(column, "") for column in columns)) + " |"
        )
    lines.extend(
        [
            "",
            "Notes:",
            "",
            "- Intervals are percentile bootstrap confidence intervals.",
            "- Repeated-CV rows use isolate-cluster bootstrap, keeping each isolate's repeated predictions together.",
            "- UNKNOWN/no-call coverage is reported separately from binary performance metrics.",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> int:
    args = parse_args()
    if args.iterations < 100:
        raise ValueError("--iterations should be at least 100 for stable intervals")
    table_rows = read_csv(args.table)
    ci_rows, table_with_ci = build_ci_rows(
        table_rows,
        read_csv(args.isolates),
        read_csv(args.five_level),
        read_csv(args.cv_predictions),
        read_csv(args.site_predictions),
        args.iterations,
        args.confidence,
        args.seed,
    )
    table_ci_fieldnames = list(table_rows[0].keys()) + [
        f"{metric}_ci" for metric in METRICS
    ]
    write_csv(args.ci_output, ci_rows, CI_FIELDNAMES)
    write_csv(args.table_ci_output, table_with_ci, table_ci_fieldnames)
    args.table_ci_md_output.parent.mkdir(parents=True, exist_ok=True)
    args.table_ci_md_output.write_text(markdown_table(table_with_ci), encoding="utf-8")
    print(f"ci={args.ci_output}")
    print(f"table_ci={args.table_ci_output}")
    print(f"table_ci_markdown={args.table_ci_md_output}")
    print(f"rows={len(table_with_ci)}")
    print(f"iterations={args.iterations}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
