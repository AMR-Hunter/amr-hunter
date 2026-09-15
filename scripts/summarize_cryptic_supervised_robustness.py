from __future__ import annotations
import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path

DEFAULT_PREDICTIONS = Path("reports/publication/cryptic_supervised_cv_predictions.csv")
DEFAULT_REPEAT_OUTPUT = Path(
    "reports/publication/cryptic_supervised_repeat_metrics.csv"
)
DEFAULT_SUMMARY_OUTPUT = Path(
    "reports/publication/cryptic_supervised_robustness_summary.csv"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize repeat-level stability for supervised CRyPTIC BDQ CV."
    )
    parser.add_argument("--predictions", type=Path, default=DEFAULT_PREDICTIONS)
    parser.add_argument("--repeat-output", type=Path, default=DEFAULT_REPEAT_OUTPUT)
    parser.add_argument("--summary-output", type=Path, default=DEFAULT_SUMMARY_OUTPUT)
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


def ratio(num: int, den: int) -> float | None:
    return num / den if den else None


def format_float(value: float | None) -> str:
    return f"{value:.4f}" if value is not None else ""


def confusion_counts(truth: str, prediction: str) -> str:
    is_r = truth == "R"
    pred_r = prediction == "R"
    if is_r and pred_r:
        return "tp"
    if is_r and (not pred_r):
        return "fn"
    if not is_r and pred_r:
        return "fp"
    return "tn"


def metric_row(
    model: str, repeat: str, target: str, counts: Counter[str]
) -> dict[str, str]:
    total = counts["tp"] + counts["tn"] + counts["fp"] + counts["fn"]
    rec = ratio(counts["tp"], counts["tp"] + counts["fn"])
    spec = ratio(counts["tn"], counts["tn"] + counts["fp"])
    prec = ratio(counts["tp"], counts["tp"] + counts["fp"])
    acc = ratio(counts["tp"] + counts["tn"], total)
    return {
        "model": model,
        "repeat": repeat,
        "target_specificity": target,
        "total": str(total),
        "tp": str(counts["tp"]),
        "tn": str(counts["tn"]),
        "fp": str(counts["fp"]),
        "fn": str(counts["fn"]),
        "recall": format_float(rec),
        "specificity": format_float(spec),
        "precision": format_float(prec),
        "accuracy": format_float(acc),
    }


def mean(values: list[float]) -> float | None:
    return sum(values) / len(values) if values else None


def sample_sd(values: list[float]) -> float | None:
    if len(values) < 2:
        return None
    avg = sum(values) / len(values)
    return math.sqrt(sum(((value - avg) ** 2 for value in values)) / (len(values) - 1))


def summarize_metric(values: list[float]) -> dict[str, str]:
    if not values:
        return {"mean": "", "sd": "", "min": "", "max": ""}
    return {
        "mean": format_float(mean(values)),
        "sd": format_float(sample_sd(values)),
        "min": format_float(min(values)),
        "max": format_float(max(values)),
    }


def target_columns(fieldnames: list[str]) -> list[tuple[str, str]]:
    columns = []
    for name in fieldnames:
        if name.startswith("pred_at_spec_"):
            columns.append((name, name.removeprefix("pred_at_spec_")))
    return columns


def build_repeat_metrics(
    rows: list[dict[str, str]], fieldnames: list[str]
) -> list[dict[str, str]]:
    groups: dict[tuple[str, str, str], Counter[str]] = defaultdict(Counter)
    for row in rows:
        for column, target in target_columns(fieldnames):
            key = (row["model"], row["repeat"], target)
            groups[key][confusion_counts(row["bdq_binary_phenotype"], row[column])] += 1
    out = []
    for model, repeat, target in sorted(groups):
        out.append(metric_row(model, repeat, target, groups[model, repeat, target]))
    return out


def build_summary(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["model"], row["target_specificity"]].append(row)
    out = []
    for model, target in sorted(grouped):
        repeat_rows = grouped[model, target]
        metrics = {
            name: [float(row[name]) for row in repeat_rows if row[name]]
            for name in ("recall", "specificity", "precision", "accuracy")
        }
        row_out = {
            "model": model,
            "target_specificity": target,
            "repeats": str(len(repeat_rows)),
        }
        for metric_name, values in metrics.items():
            summary = summarize_metric(values)
            for stat, value in summary.items():
                row_out[f"{metric_name}_{stat}"] = value
        out.append(row_out)
    return out


def main() -> int:
    args = parse_args()
    rows, fieldnames = read_csv(args.predictions)
    repeat_rows = build_repeat_metrics(rows, fieldnames)
    summary_rows = build_summary(repeat_rows)
    write_csv(args.repeat_output, repeat_rows)
    write_csv(args.summary_output, summary_rows)
    print(f"repeat_metrics={args.repeat_output}")
    print(f"repeat_rows={len(repeat_rows)}")
    print(f"summary={args.summary_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
