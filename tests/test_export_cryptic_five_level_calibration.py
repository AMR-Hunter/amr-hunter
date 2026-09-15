from collections import Counter
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
from scripts.export_cryptic_five_level_calibration import (
    assign_five_level_label,
    build_carrier_stats,
    build_label_summary,
    build_metric_rows,
)


def isolate(unique_id: str, phenotype: str, prediction: str) -> dict[str, str]:
    return {
        "dataset": "CRyPTIC",
        "unique_id": unique_id,
        "ena_sample": unique_id,
        "bdq_binary_phenotype": phenotype,
        "bdq_mic": "0.5" if phenotype == "R" else "0.06",
        "amr_hunter_isolate_prediction": prediction,
        "call_source_summary": "fixture",
    }


def variant(
    unique_id: str,
    phenotype: str,
    gene: str,
    mutation: str,
    prediction: str,
    source: str,
) -> dict[str, str]:
    return {
        "unique_id": unique_id,
        "bdq_binary_phenotype": phenotype,
        "gene": gene,
        "mutation": mutation,
        "amr_hunter_prediction": prediction,
        "call_source": source,
    }


def test_five_level_labels_and_metrics():
    isolates = [
        isolate("strong-r", "R", "R"),
        isolate("likely-r", "R", "R"),
        isolate("weak-r", "S", "R"),
        isolate("core-unknown", "S", "S"),
        isolate("context", "S", "S"),
        isolate("plain-s", "S", "S"),
        isolate("unknown", "S", "UNKNOWN"),
    ]
    variants = [
        variant("strong-r", "R", "Rv0678", "p.Gly155fs", "R", "consequence_inferred"),
        variant(
            "strong-r-sensitive-carrier",
            "S",
            "Rv0678",
            "p.Gly155fs",
            "R",
            "consequence_inferred",
        ),
        variant("likely-r", "R", "Rv0678", "p.Asp47fs", "R", "consequence_inferred"),
        variant("weak-r", "S", "glpK", "p.Val353Phe", "R", "snapshot_virtual"),
        variant(
            "core-unknown",
            "S",
            "Rv0678",
            "p.Arg90Cys",
            "UNKNOWN",
            "unsupported_variant",
        ),
        variant(
            "context", "S", "mtrB", "p.Met517Leu", "UNKNOWN", "unsupported_variant"
        ),
        variant("plain-s", "S", "mmpL5", "p.Ile948Val", "S", "direct_record_match"),
        variant("carrier-a", "S", "Rv0678", "p.Asp47fs", "R", "consequence_inferred"),
        variant("carrier-b", "S", "Rv0678", "p.Asp47fs", "R", "consequence_inferred"),
        variant("carrier-c", "S", "Rv0678", "p.Asp47fs", "R", "consequence_inferred"),
        variant("carrier-d", "S", "Rv0678", "p.Asp47fs", "R", "consequence_inferred"),
    ]
    carrier_stats = build_carrier_stats(variants)
    by_isolate = {}
    for row in variants:
        by_isolate.setdefault(row["unique_id"], []).append(row)
    labels = [
        assign_five_level_label(
            row, by_isolate.get(row["unique_id"], []), carrier_stats
        )
        for row in isolates
    ]
    by_id = {row["unique_id"]: row for row in labels}
    assert by_id["strong-r"]["five_level_label"] == "R"
    assert by_id["likely-r"]["five_level_label"] == "Likely R"
    assert by_id["weak-r"]["five_level_label"] == "Indeterminate"
    assert by_id["core-unknown"]["five_level_label"] == "Indeterminate"
    assert by_id["context"]["five_level_label"] == "Likely S"
    assert by_id["plain-s"]["five_level_label"] == "S"
    assert by_id["unknown"]["five_level_label"] == "Indeterminate"
    summary = build_label_summary(labels)
    summary_counts = {
        row["label"]: (row["phenotype_r"], row["phenotype_s"]) for row in summary
    }
    assert summary_counts["R"] == ("1", "0")
    assert summary_counts["Likely R"] == ("1", "0")
    assert summary_counts["Indeterminate"] == ("0", "3")
    metrics = {row["metric"]: row for row in build_metric_rows(labels)}
    assert metrics["hard_r_precision"]["value"] == "1.0000"
    assert metrics["r_or_likely_r_recall"]["value"] == "1.0000"
    assert metrics["indeterminate_rate"]["value"] == "0.4286"


if __name__ == "__main__":
    test_five_level_labels_and_metrics()
    print("cryptic five-level calibration tests passed")
