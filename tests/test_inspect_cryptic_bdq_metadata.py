import csv
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURE = REPO_ROOT / "tests" / "fixtures" / "cryptic_reuse_table_minimal.csv"
SCRIPT = REPO_ROOT / "scripts" / "inspect_cryptic_bdq_metadata.py"


def test_inspection_script_writes_expected_metrics(tmp_path: Path) -> None:
    output = tmp_path / "metadata_inspection.csv"
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--input", str(FIXTURE), "--output", str(output)],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=True,
    )
    assert "validation_feasibility=feasible" in result.stdout
    with output.open(newline="", encoding="utf-8") as handle:
        rows = {row["metric"]: row["value"] for row in csv.DictReader(handle)}
    assert rows["total_rows"] == "5"
    assert rows["bdq_binary_phenotype_r"] == "2"
    assert rows["bdq_binary_phenotype_s"] == "2"
    assert rows["bdq_binary_phenotype_missing"] == "1"
    assert rows["bdq_mic_available"] == "4"
    assert rows["vcf_available"] == "4"
    assert rows["quality_passing_r"] == "1"
    assert rows["quality_passing_s"] == "2"
    assert rows["quality_passing_with_vcf_r"] == "1"
    assert rows["quality_passing_with_vcf_s"] == "1"


if __name__ == "__main__":
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        test_inspection_script_writes_expected_metrics(Path(tmp))
    print("cryptic metadata inspection tests passed")
