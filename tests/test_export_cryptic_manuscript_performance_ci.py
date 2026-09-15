from pathlib import Path
import random
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
from scripts.export_cryptic_manuscript_performance_ci import (
    bootstrap_metric_cis,
    counts_for,
    five_level_groups,
    frozen_groups,
    metric_values,
    prediction_groups,
)


def test_counts_and_metric_values():
    assert counts_for("R", "R") == (1, 0, 0, 0)
    assert counts_for("S", "S") == (0, 1, 0, 0)
    assert counts_for("S", "R") == (0, 0, 1, 0)
    assert counts_for("R", "S") == (0, 0, 0, 1)
    assert counts_for("S", None) == (0, 0, 0, 0)
    values = metric_values((2, 6, 2, 1))
    assert round(values["recall"], 4) == 0.6667
    assert values["specificity"] == 0.75
    assert values["precision"] == 0.5


def test_group_builders_match_binary_policies():
    isolates = [
        {"bdq_binary_phenotype": "R", "amr_hunter_isolate_prediction": "R"},
        {"bdq_binary_phenotype": "S", "amr_hunter_isolate_prediction": "UNKNOWN"},
    ]
    assert frozen_groups(isolates) == [(1, 0, 0, 0), (0, 0, 0, 0)]
    five = [
        {"bdq_binary_phenotype": "R", "five_level_label": "Likely R"},
        {"bdq_binary_phenotype": "S", "five_level_label": "Indeterminate"},
    ]
    assert five_level_groups(five, {"R"}) == [(0, 0, 0, 1), (0, 1, 0, 0)]
    assert five_level_groups(five, {"R", "Likely R"}) == [(1, 0, 0, 0), (0, 1, 0, 0)]


def test_prediction_groups_can_cluster_by_isolate():
    rows = [
        {
            "model": "m",
            "unique_id": "iso1",
            "bdq_binary_phenotype": "R",
            "pred_at_spec_0.95": "R",
        },
        {
            "model": "m",
            "unique_id": "iso1",
            "bdq_binary_phenotype": "R",
            "pred_at_spec_0.95": "S",
        },
        {
            "model": "m",
            "unique_id": "iso2",
            "bdq_binary_phenotype": "S",
            "pred_at_spec_0.95": "S",
        },
    ]
    clustered = prediction_groups(rows, "m", "0.95", cluster_by_isolate=True)
    unclustered = prediction_groups(rows, "m", "0.95", cluster_by_isolate=False)
    assert sorted(clustered) == sorted([(1, 0, 0, 1), (0, 1, 0, 0)])
    assert unclustered == [(1, 0, 0, 0), (0, 0, 0, 1), (0, 1, 0, 0)]


def test_bootstrap_metric_cis_returns_bounds_for_available_metrics():
    groups = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0)]
    cis = bootstrap_metric_cis(
        groups, iterations=100, confidence=0.95, rng=random.Random(1)
    )
    assert 0 <= cis["specificity"][0] <= cis["specificity"][1] <= 1
    assert 0 <= cis["precision"][0] <= cis["precision"][1] <= 1
    assert cis["recall"][0] == 1.0
    assert cis["recall"][1] == 1.0


if __name__ == "__main__":
    test_counts_and_metric_values()
    test_group_builders_match_binary_policies()
    test_prediction_groups_can_cluster_by_isolate()
    test_bootstrap_metric_cis_returns_bounds_for_available_metrics()
    print("CRyPTIC manuscript performance CI tests passed")
