import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
from scripts.cryptic_bdq_validation_common import (
    ValidationMetrics,
    format_ratio,
    normalize_bdq_phenotype,
    passes_bdq_quality,
    summarize_metadata_rows,
)


def test_bdq_phenotype_normalization() -> None:
    assert normalize_bdq_phenotype("R") == "R"
    assert normalize_bdq_phenotype("resistant") == "R"
    assert normalize_bdq_phenotype("susceptible") == "S"
    assert normalize_bdq_phenotype("mixed") == ""


def test_quality_filter_marks_failed_and_missing_as_not_passing() -> None:
    assert passes_bdq_quality("HIGH")
    assert passes_bdq_quality("pass")
    assert not passes_bdq_quality("failed")
    assert not passes_bdq_quality("")


def test_ratio_and_validation_metrics_denominators() -> None:
    assert format_ratio(1, 4) == "0.2500"
    assert format_ratio(1, 0) == ""
    metrics = ValidationMetrics(
        total_isolates=5,
        phenotype_r=2,
        phenotype_s=3,
        interpretable_isolates=4,
        unknown_isolates=1,
        tp=1,
        tn=2,
        fp=1,
        fn=0,
    )
    assert metrics.coverage == "0.8000"
    assert metrics.accuracy == "0.7500"
    assert metrics.sensitivity == "1.0000"
    assert metrics.specificity == "0.6667"


def test_metadata_summary_counts_quality_and_vcf_denominators() -> None:
    rows = [
        {
            "BDQ_BINARY_PHENOTYPE": "R",
            "BDQ_MIC": "1.0",
            "BDQ_PHENOTYPE_QUALITY": "high",
            "VCF": "a.vcf.gz",
        },
        {
            "BDQ_BINARY_PHENOTYPE": "S",
            "BDQ_MIC": "0.03",
            "BDQ_PHENOTYPE_QUALITY": "pass",
            "VCF": "b.vcf.gz",
        },
        {
            "BDQ_BINARY_PHENOTYPE": "R",
            "BDQ_MIC": "",
            "BDQ_PHENOTYPE_QUALITY": "failed",
            "VCF": "c.vcf.gz",
        },
    ]
    summary = {row["metric"]: row["value"] for row in summarize_metadata_rows(rows)}
    assert summary["total_rows"] == "3"
    assert summary["bdq_binary_phenotype_r"] == "2"
    assert summary["bdq_binary_phenotype_s"] == "1"
    assert summary["bdq_mic_available"] == "2"
    assert summary["quality_passing_with_vcf_r"] == "1"
    assert summary["quality_passing_with_vcf_s"] == "1"
    assert summary["high_quality_with_vcf_r"] == "1"
    assert summary["high_quality_with_vcf_s"] == "0"
    assert summary["validation_feasibility"] == "feasible"


if __name__ == "__main__":
    test_bdq_phenotype_normalization()
    test_quality_filter_marks_failed_and_missing_as_not_passing()
    test_ratio_and_validation_metrics_denominators()
    test_metadata_summary_counts_quality_and_vcf_denominators()
    print("cryptic validation common tests passed")
