from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
from scripts.run_cryptic_supervised_site_holdout import (
    add_site,
    feature_names,
    infer_site,
    run_site_holdout,
)


def row(
    unique_id: str, phenotype: str, signal: str, baseline: str = "S"
) -> dict[str, str]:
    return {
        "unique_id": unique_id,
        "ena_sample": unique_id,
        "bdq_binary_phenotype": phenotype,
        "label": "1" if phenotype == "R" else "0",
        "bdq_mic": "0.5" if phenotype == "R" else "0.06",
        "baseline_prediction": baseline,
        "baseline_pred_R": "1" if baseline == "R" else "0",
        "baseline_pred_S": "1" if baseline == "S" else "0",
        "signal": signal,
    }


def test_infer_site_and_feature_names():
    assert infer_site("site.02.subj.0001.lab.X.iso.1") == "site.02"
    assert infer_site("missing-prefix") == "unknown"
    fields = [
        "unique_id",
        "label",
        "bdq_mic",
        "baseline_prediction",
        "site",
        "baseline_pred_R",
        "signal",
    ]
    assert feature_names(fields) == ["baseline_pred_R", "signal"]


def test_run_site_holdout_returns_group_and_summary_metrics():
    rows = add_site(
        [
            row("site.01.iso.r1", "R", "1", baseline="R"),
            row("site.01.iso.s1", "S", "0"),
            row("site.02.iso.r1", "R", "1", baseline="R"),
            row("site.02.iso.s1", "S", "0"),
            row("site.03.iso.r1", "R", "1", baseline="R"),
            row("site.03.iso.s1", "S", "0"),
        ]
    )
    summary_rows, metric_rows, prediction_rows = run_site_holdout(
        rows, ["baseline_pred_R", "baseline_pred_S", "signal"]
    )
    summary_by_key = {
        (row["model"], row["target_specificity"]): row for row in summary_rows
    }
    baseline = summary_by_key["baseline_only", "0.98"]
    assert baseline["heldout_groups"] == "3"
    assert baseline["recall"] == "1.0000"
    assert baseline["specificity"] == "1.0000"
    assert {row["heldout_site"] for row in metric_rows} == {
        "site.01",
        "site.02",
        "site.03",
    }
    assert len(prediction_rows) == 6 * 5


if __name__ == "__main__":
    test_infer_site_and_feature_names()
    test_run_site_holdout_returns_group_and_summary_metrics()
    print("CRyPTIC supervised site holdout tests passed")
