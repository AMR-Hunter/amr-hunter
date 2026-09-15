import csv
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "extract_cryptic_bdq_scope_variants.py"
METADATA = REPO_ROOT / "tests" / "fixtures" / "cryptic_adapter_metadata.csv"
REFERENCE = REPO_ROOT / "tests" / "fixtures" / "cryptic_adapter_reference.gbk"
VCF_ROOT = REPO_ROOT / "tests" / "fixtures"


def test_extract_adapter_writes_normalized_variant_rows(tmp_path: Path) -> None:
    output = tmp_path / "bdq_scope_variant_calls.tsv"
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--metadata",
            str(METADATA),
            "--reference-gbk",
            str(REFERENCE),
            "--vcf-root",
            str(VCF_ROOT),
            "--output",
            str(output),
            "--genes",
            "atpE,Rv0678",
            "--max-isolates",
            "1",
        ],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=True,
    )
    assert "eligible_isolates=1" in result.stdout
    with output.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["gene"] == "atpE"
    assert rows[0]["mutation"] == "p.Ala2Pro"
    assert rows[0]["region_type"] == "CDS"
    assert rows[0]["consequence"] == "missense_variant"
    assert rows[1]["gene"] == "atpE"
    assert rows[1]["mutation"] == "c.-6A>G"
    assert rows[1]["region_type"] == "REGULATORY"
    assert rows[2]["variant_class"] == "frameshift"
    assert rows[2]["consequence"] == "frameshift_variant"
    assert len(rows) == 3


if __name__ == "__main__":
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        test_extract_adapter_writes_normalized_variant_rows(Path(tmp))
    print("cryptic VCF adapter tests passed")
