from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
from scripts.summarize_cryptic_supervised_robustness import (
    build_repeat_metrics,
    build_summary,
    target_columns,
)


def test_build_repeat_metrics_and_summary():
    fieldnames = [
        "model",
        "repeat",
        "bdq_binary_phenotype",
        "pred_at_spec_0.95",
        "pred_at_spec_0.90",
    ]
    rows = [
        {
            "model": "all_features",
            "repeat": "1",
            "bdq_binary_phenotype": "R",
            "pred_at_spec_0.95": "R",
            "pred_at_spec_0.90": "R",
        },
        {
            "model": "all_features",
            "repeat": "1",
            "bdq_binary_phenotype": "S",
            "pred_at_spec_0.95": "S",
            "pred_at_spec_0.90": "R",
        },
        {
            "model": "all_features",
            "repeat": "2",
            "bdq_binary_phenotype": "R",
            "pred_at_spec_0.95": "S",
            "pred_at_spec_0.90": "R",
        },
        {
            "model": "all_features",
            "repeat": "2",
            "bdq_binary_phenotype": "S",
            "pred_at_spec_0.95": "S",
            "pred_at_spec_0.90": "S",
        },
    ]
    assert target_columns(fieldnames) == [
        ("pred_at_spec_0.95", "0.95"),
        ("pred_at_spec_0.90", "0.90"),
    ]
    repeat_rows = build_repeat_metrics(rows, fieldnames)
    by_key = {(row["repeat"], row["target_specificity"]): row for row in repeat_rows}
    assert by_key["1", "0.95"]["tp"] == "1"
    assert by_key["1", "0.95"]["specificity"] == "1.0000"
    assert by_key["1", "0.90"]["fp"] == "1"
    assert by_key["2", "0.95"]["recall"] == "0.0000"
    summary = {row["target_specificity"]: row for row in build_summary(repeat_rows)}
    assert summary["0.95"]["repeats"] == "2"
    assert summary["0.95"]["recall_mean"] == "0.5000"
    assert summary["0.90"]["recall_mean"] == "1.0000"


if __name__ == "__main__":
    test_build_repeat_metrics_and_summary()
    print("CRyPTIC supervised robustness tests passed")
