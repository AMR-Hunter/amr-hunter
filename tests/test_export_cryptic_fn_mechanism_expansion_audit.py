from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
from scripts.export_cryptic_fn_mechanism_expansion_audit import (
    build_fn_mic_rows,
    build_fn_mic_summary_rows,
    build_indices,
    combination_enrichment_rows,
    combo_memberships,
    fisher_one_sided_enrichment,
    odds_ratio_haldane,
    priority_rows,
    recurrent_fn_variant_rows,
    variant_enrichment_rows,
)


def isolate(
    unique_id: str, phenotype: str, prediction: str, mic: str = "0.5"
) -> dict[str, str]:
    return {
        "unique_id": unique_id,
        "ena_sample": unique_id,
        "bdq_binary_phenotype": phenotype,
        "bdq_mic": mic,
        "amr_hunter_isolate_prediction": prediction,
    }


def variant(
    unique_id: str,
    gene: str,
    mutation: str,
    source: str = "unsupported_variant",
    prediction: str = "UNKNOWN",
) -> dict[str, str]:
    return {
        "unique_id": unique_id,
        "gene": gene,
        "mutation": mutation,
        "call_source": source,
        "amr_hunter_prediction": prediction,
    }


def test_enrichment_helpers_are_directional():
    enriched_p = fisher_one_sided_enrichment(a=3, b=0, c=1, d=8)
    weak_p = fisher_one_sided_enrichment(a=1, b=5, c=3, d=3)
    assert enriched_p < weak_p
    assert odds_ratio_haldane(3, 0, 1, 8) > odds_ratio_haldane(1, 5, 3, 3)


def test_fn_rows_enrichment_and_priority_candidates():
    isolates = [
        isolate("r-fn-1", "R", "S", "0.5"),
        isolate("r-fn-2", "R", "S", "1.0"),
        isolate("r-tp", "R", "R", "2.0"),
        isolate("s-1", "S", "S", "0.06"),
        isolate("s-2", "S", "S", "0.06"),
        isolate("s-3", "S", "S", "0.06"),
    ]
    variants = [
        variant("r-fn-1", "Rv0678", "p.Ala10Val"),
        variant("r-fn-1", "mtrB", "p.Met517Leu"),
        variant("r-fn-1", "glpK", "p.Val192Met"),
        variant("r-fn-2", "Rv0678", "p.Ala10Val"),
        variant("r-fn-2", "mtrB", "p.Met517Leu"),
        variant("r-tp", "Rv0678", "p.Ala10Val"),
        variant("s-1", "mtrB", "p.Met517Leu"),
        variant("s-2", "mtrB", "p.Met517Leu"),
    ]
    _isolates_by_id, variants_by_id = build_indices(isolates, variants)
    fn_rows = build_fn_mic_rows(isolates, variants_by_id)
    fn_summary = build_fn_mic_summary_rows(fn_rows)
    assert len(fn_rows) == 2
    assert fn_summary[0]["fn_isolates"] == "1"
    assert fn_rows[0]["unsupported_gene_count"] == "3"
    assert "Rv0678:p.Ala10Val" in fn_rows[0]["unsupported_variants"]
    memberships = combo_memberships(variants_by_id["r-fn-1"])
    assert "gene_pair:Rv0678+glpK" in memberships
    assert "context:mtrB_Met517Leu_plus_secondary_gene" in memberships
    assert "context:mtrB_Met517Leu_plus_core_gene" in memberships
    enrichment = variant_enrichment_rows(isolates, variants_by_id, min_carriers=1)
    by_variant = {row["variant"]: row for row in enrichment}
    assert by_variant["Rv0678:p.Ala10Val"]["r_carriers"] == "3"
    assert by_variant["Rv0678:p.Ala10Val"]["s_carriers"] == "0"
    assert by_variant["mtrB:p.Met517Leu"]["is_common_background"] == "true"
    recurrent = recurrent_fn_variant_rows(enrichment, min_fn_carriers=2)
    assert {row["variant"] for row in recurrent} == {
        "Rv0678:p.Ala10Val",
        "mtrB:p.Met517Leu",
    }
    combos = combination_enrichment_rows(isolates, variants_by_id, min_carriers=1)
    priority = priority_rows(recurrent, combos)
    assert priority[0]["priority_tier"] == "high_priority_candidate"
    assert any((row["candidate"] == "mtrB:p.Met517Leu" for row in priority))


if __name__ == "__main__":
    test_enrichment_helpers_are_directional()
    test_fn_rows_enrichment_and_priority_candidates()
    print("CRyPTIC FN mechanism expansion audit tests passed")
