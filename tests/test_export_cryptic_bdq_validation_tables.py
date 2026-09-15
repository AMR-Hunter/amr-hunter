import csv
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "export_cryptic_bdq_validation_tables.py"
METADATA = REPO_ROOT / "tests" / "fixtures" / "cryptic_validation_metadata_minimal.csv"
VARIANTS = REPO_ROOT / "tests" / "fixtures" / "cryptic_variant_calls_minimal.tsv"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def test_exporter_writes_summary_isolate_and_variant_tables(tmp_path: Path) -> None:
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--metadata",
            str(METADATA),
            "--variant-calls",
            str(VARIANTS),
            "--quality-filter",
            "non_missing",
            "--output-dir",
            str(tmp_path),
        ],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=True,
    )
    assert "total=6 interpretable=4 coverage=0.6667 accuracy=1.0000" in result.stdout
    summary = read_csv(tmp_path / "external_bdq_phenotype_validation_summary.csv")
    isolates = read_csv(tmp_path / "external_bdq_phenotype_validation_isolates.csv")
    variants = read_csv(
        tmp_path / "external_bdq_phenotype_validation_variant_calls.csv"
    )
    assert summary == [
        {
            "dataset": "CRyPTIC",
            "filter_name": "non_missing_quality_with_vcf",
            "total_isolates": "6",
            "phenotype_r": "2",
            "phenotype_s": "4",
            "interpretable_isolates": "4",
            "unknown_isolates": "2",
            "coverage": "0.6667",
            "accuracy": "1.0000",
            "tp": "2",
            "tn": "2",
            "fp": "0",
            "fn": "0",
            "sensitivity": "1.0000",
            "specificity": "1.0000",
            "notes": "Accuracy excludes UNKNOWN/no-call isolates; coverage includes all included isolates.",
        }
    ]
    by_isolate = {row["unique_id"]: row for row in isolates}
    assert by_isolate["iso-r-direct"]["amr_hunter_isolate_prediction"] == "R"
    assert by_isolate["iso-r-direct"]["top_resistance_variant"] == "atpE:p.Ala63Pro"
    assert by_isolate["iso-r-lof"]["call_source_summary"] == "consequence_inferred:1"
    assert by_isolate["iso-s-snapshot"]["call_source_summary"] == "snapshot_virtual:1"
    assert by_isolate["iso-s-unsupported"]["amr_hunter_isolate_prediction"] == "UNKNOWN"
    assert (
        by_isolate["iso-s-unsupported"]["unknown_reason"]
        == "unsupported_bdq_scope_variant"
    )
    assert by_isolate["iso-s-nocall"]["unknown_reason"] == "no_bdq_scope_variant"
    call_sources = {row["unique_id"]: row["call_source"] for row in variants}
    assert call_sources["iso-r-direct"] == "direct_record_match"
    assert call_sources["iso-r-lof"] == "consequence_inferred"
    assert call_sources["iso-s-snapshot"] == "snapshot_virtual"
    assert call_sources["iso-s-unsupported"] == "unsupported_variant"
    assert call_sources["iso-s-nocall"] == "no_call"


if __name__ == "__main__":
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        test_exporter_writes_summary_isolate_and_variant_tables(Path(tmp))
    print("cryptic validation exporter tests passed")
