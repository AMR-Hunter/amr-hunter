from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
from scripts.export_cryptic_mtrb_enrichment import (
    build_combo_summary,
    build_group_summary,
    group_memberships,
)


def isolate(unique_id: str, phenotype: str, prediction: str = "S") -> dict[str, str]:
    return {
        "unique_id": unique_id,
        "ena_sample": unique_id,
        "bdq_binary_phenotype": phenotype,
        "bdq_mic": "0.5" if phenotype == "R" else "0.06",
        "amr_hunter_isolate_prediction": prediction,
        "call_source_summary": "fixture",
    }


def variant(
    unique_id: str,
    gene: str,
    mutation: str,
    prediction: str = "UNKNOWN",
    source: str = "unsupported_variant",
) -> dict[str, str]:
    return {
        "unique_id": unique_id,
        "gene": gene,
        "mutation": mutation,
        "amr_hunter_prediction": prediction,
        "call_source": source,
    }


def test_group_memberships_and_enrichment_summaries():
    rows = [
        variant("iso-r", "mtrB", "p.Met517Leu"),
        variant("iso-r", "mtrB", "p.Pro18Ser"),
        variant("iso-r", "Rv0678", "p.Arg90Cys"),
        variant("iso-r", "glpK", "p.Val192Met"),
    ]
    memberships = group_memberships(rows)
    assert "mtrB_context_any" in memberships
    assert "mtrB_double_hit" in memberships
    assert "mtrB_Met517Leu_plus_core_unsupported" in memberships
    assert "mtrB_double_hit_plus_secondary_unsupported" in memberships
    isolates = [
        isolate("iso-r", "R"),
        isolate("iso-s", "S"),
        isolate("iso-no-mtrb", "S"),
    ]
    variants_by_isolate = {
        "iso-r": rows,
        "iso-s": [
            variant("iso-s", "mtrB", "p.Met517Leu"),
            variant("iso-s", "mtrB", "p.Pro18Ser"),
        ],
        "iso-no-mtrb": [variant("iso-no-mtrb", "Rv0678", "p.Arg90Cys")],
    }
    five_level = {
        "iso-r": {"five_level_label": "Indeterminate"},
        "iso-s": {"five_level_label": "Likely S"},
        "iso-no-mtrb": {"five_level_label": "Indeterminate"},
    }
    groups = {
        row["group"]: row
        for row in build_group_summary(isolates, variants_by_isolate, five_level)
    }
    assert groups["mtrB_double_hit"]["total_isolates"] == "2"
    assert groups["mtrB_double_hit"]["phenotype_r"] == "1"
    assert groups["mtrB_context_absent"]["total_isolates"] == "1"
    combos = build_combo_summary(isolates, variants_by_isolate, five_level)
    top = combos[0]
    assert "p.Met517Leu" in top["mtrB_signature"]
    assert "p.Pro18Ser" in top["mtrB_signature"]
    assert top["phenotype_r"] == "1"
    assert "glpK:p.Val192Met:1" in top["top_secondary_unsupported"]
    assert "Rv0678:p.Arg90Cys:1" in top["top_core_unsupported"]


if __name__ == "__main__":
    test_group_memberships_and_enrichment_summaries()
    print("cryptic mtrB enrichment tests passed")
