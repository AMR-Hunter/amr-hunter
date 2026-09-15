from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
from scripts.export_cryptic_candidate_nested_validation import (
    NESTED_MODEL,
    REFERENCE_MODEL,
    unsupported_variants_by_isolate,
)
from scripts.run_cryptic_candidate_nested_site_holdout import (
    add_site,
    infer_site,
    run_site_holdout_nested,
)


def row(unique_id: str, phenotype: str, baseline: str = "S") -> dict[str, str]:
    return {
        "unique_id": unique_id,
        "ena_sample": unique_id,
        "bdq_binary_phenotype": phenotype,
        "label": "1" if phenotype == "R" else "0",
        "bdq_mic": "0.5" if phenotype == "R" else "0.06",
        "baseline_prediction": baseline,
        "baseline_pred_R": "1" if baseline == "R" else "0",
        "baseline_pred_S": "1" if baseline == "S" else "0",
        "secondary_unsupported": "0",
    }


def variant(unique_id: str, gene: str, mutation: str) -> dict[str, str]:
    return {
        "unique_id": unique_id,
        "gene": gene,
        "mutation": mutation,
        "call_source": "unsupported_variant",
    }


def test_infer_site():
    assert infer_site("site.02.subj.0001.lab.x.iso.1") == "site.02"
    assert infer_site("no-site") == "unknown"


def test_run_site_holdout_nested_outputs_expected_rows():
    rows = add_site(
        [
            row("site.01.r1", "R"),
            row("site.01.s1", "S"),
            row("site.02.r1", "R"),
            row("site.02.s1", "S"),
            row("site.03.r1", "R"),
            row("site.03.s1", "S"),
        ]
    )
    variants = [
        variant("site.01.r1", "Rv1979c", "p.Arg409Gln"),
        variant("site.02.r1", "Rv1979c", "p.Arg409Gln"),
        variant("site.03.r1", "Rv1979c", "p.Arg409Gln"),
        variant("site.01.s1", "mtrB", "p.Met517Leu"),
    ]
    unsupported = unsupported_variants_by_isolate(variants)
    summary, metrics, predictions, selections = run_site_holdout_nested(
        rows,
        ["baseline_pred_R", "baseline_pred_S", "secondary_unsupported"],
        unsupported,
        min_train_r_carriers=1,
        min_train_fn_carriers=1,
        min_r_fraction=0.05,
        min_odds_ratio=1.0,
        max_p_value=1.0,
        min_train_positive=1,
    )
    assert len(summary) == 8
    assert len(metrics) == 3 * 2 * 4
    assert len(predictions) == len(rows) * 2
    assert selections
    assert {row["model"] for row in summary} == {REFERENCE_MODEL, NESTED_MODEL}
    assert {row["heldout_site"] for row in metrics} == {"site.01", "site.02", "site.03"}


if __name__ == "__main__":
    test_infer_site()
    test_run_site_holdout_nested_outputs_expected_rows()
    print("CRyPTIC candidate nested site holdout tests passed")
