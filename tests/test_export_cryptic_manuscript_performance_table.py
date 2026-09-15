from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
from scripts.export_cryptic_manuscript_performance_table import (
    binary_row,
    build_rows,
    five_level_row,
    markdown_table,
    selected_supervised_row,
)


def test_binary_row_adds_precision_and_preserves_coverage():
    rows = [
        {
            "total_isolates": "10",
            "phenotype_r": "4",
            "phenotype_s": "6",
            "interpretable_isolates": "9",
            "unknown_isolates": "1",
            "coverage": "0.9000",
            "accuracy": "0.7778",
            "tp": "2",
            "tn": "5",
            "fp": "1",
            "fn": "1",
            "sensitivity": "0.6667",
            "specificity": "0.8333",
        }
    ]
    out = binary_row(rows)
    assert out["coverage"] == "0.9000"
    assert out["precision"] == "0.6667"
    assert out["unknown_no_call"] == "1"


def test_five_level_projection_keeps_indeterminate_visible():
    label_rows = [
        {"label": "R", "total_isolates": "2", "phenotype_r": "1", "phenotype_s": "1"},
        {
            "label": "Likely R",
            "total_isolates": "2",
            "phenotype_r": "1",
            "phenotype_s": "1",
        },
        {
            "label": "Indeterminate",
            "total_isolates": "1",
            "phenotype_r": "0",
            "phenotype_s": "1",
        },
        {
            "label": "Likely S",
            "total_isolates": "2",
            "phenotype_r": "1",
            "phenotype_s": "1",
        },
        {"label": "S", "total_isolates": "1", "phenotype_r": "0", "phenotype_s": "1"},
    ]
    hard = five_level_row(label_rows, {"R"}, "hard", "policy")
    likely = five_level_row(label_rows, {"R", "Likely R"}, "likely", "policy")
    assert hard["tp"] == "1"
    assert hard["fp"] == "1"
    assert hard["recall"] == "0.3333"
    assert hard["indeterminate"] == "1"
    assert likely["tp"] == "2"
    assert likely["fp"] == "2"
    assert likely["specificity"] == "0.6000"


def test_selected_supervised_row_and_build_rows():
    supervised = [
        {
            "model": "no_secondary",
            "target_specificity": "0.95",
            "total_predictions": "20",
            "tp": "3",
            "tn": "13",
            "fp": "2",
            "fn": "2",
            "recall": "0.6000",
            "specificity": "0.8667",
            "precision": "0.6000",
            "accuracy": "0.8000",
        }
    ]
    row = selected_supervised_row(
        supervised, "no_secondary", "0.95", "analysis", "design", "notes"
    )
    assert row["phenotype_r"] == "5"
    assert row["phenotype_s"] == "15"
    assert row["coverage"] == "1.0000"
    binary = [
        {
            "total_isolates": "10",
            "phenotype_r": "4",
            "phenotype_s": "6",
            "interpretable_isolates": "10",
            "unknown_isolates": "0",
            "coverage": "1.0000",
            "accuracy": "0.8000",
            "tp": "2",
            "tn": "6",
            "fp": "0",
            "fn": "2",
            "sensitivity": "0.5000",
            "specificity": "1.0000",
        }
    ]
    five = [
        {"label": "R", "total_isolates": "2", "phenotype_r": "1", "phenotype_s": "1"},
        {
            "label": "Likely R",
            "total_isolates": "2",
            "phenotype_r": "1",
            "phenotype_s": "1",
        },
        {
            "label": "Indeterminate",
            "total_isolates": "1",
            "phenotype_r": "0",
            "phenotype_s": "1",
        },
        {
            "label": "Likely S",
            "total_isolates": "2",
            "phenotype_r": "1",
            "phenotype_s": "1",
        },
        {"label": "S", "total_isolates": "1", "phenotype_r": "0", "phenotype_s": "1"},
    ]
    cv = supervised + [
        {**supervised[0], "model": "no_mtrb", "target_specificity": "0.90"}
    ]
    site = [
        {
            key: value
            for key, value in row.items()
            if key
            in {
                "model",
                "target_specificity",
                "tp",
                "tn",
                "fp",
                "fn",
                "recall",
                "specificity",
                "precision",
                "accuracy",
            }
        }
        for row in cv
    ]
    for row_site in site:
        row_site["total"] = "20"
    out = build_rows(binary, five, cv, site)
    assert len(out) == 7
    assert "Supervised rows" in markdown_table(out)


if __name__ == "__main__":
    test_binary_row_adds_precision_and_preserves_coverage()
    test_five_level_projection_keeps_indeterminate_visible()
    test_selected_supervised_row_and_build_rows()
    print("CRyPTIC manuscript performance table tests passed")
