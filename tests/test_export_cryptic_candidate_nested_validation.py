from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
from scripts.export_cryptic_candidate_nested_validation import (
    NESTED_MODEL,
    add_candidate_columns,
    candidate_feature_name,
    run_nested_cv,
    select_train_candidates,
    unsupported_variants_by_isolate,
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


def test_candidate_feature_name_and_columns():
    rows = [row("r1", "R"), row("s1", "S")]
    unsupported = {"r1": {"Rv1979c:p.Arg409Gln"}, "s1": set()}
    feature = candidate_feature_name("Rv1979c:p.Arg409Gln")
    assert feature == "candidate_variant__Rv1979c__p_Arg409Gln"
    added = add_candidate_columns(rows, ["Rv1979c:p.Arg409Gln"], unsupported)
    assert added == [feature]
    assert rows[0][feature] == "1"
    assert rows[1][feature] == "0"


def test_select_train_candidates_excludes_common_background():
    rows = [
        row("r1", "R"),
        row("r2", "R"),
        row("s1", "S"),
        row("s2", "S"),
        row("s3", "S"),
    ]
    variants = [
        variant("r1", "Rv1979c", "p.Arg409Gln"),
        variant("r2", "Rv1979c", "p.Arg409Gln"),
        variant("r1", "mtrB", "p.Met517Leu"),
        variant("r2", "mtrB", "p.Met517Leu"),
        variant("s1", "mtrB", "p.Met517Leu"),
        variant("s2", "mtrB", "p.Met517Leu"),
        variant("s3", "mtrB", "p.Met517Leu"),
    ]
    unsupported = unsupported_variants_by_isolate(variants)
    selected = select_train_candidates(
        rows,
        list(range(len(rows))),
        unsupported,
        min_train_r_carriers=2,
        min_train_fn_carriers=1,
        min_r_fraction=0.05,
        min_odds_ratio=2.0,
        max_p_value=0.2,
    )
    assert [candidate["candidate"] for candidate in selected] == ["Rv1979c:p.Arg409Gln"]


def test_run_nested_cv_outputs_predictions_and_selection_rows():
    rows = [
        row("r1", "R"),
        row("r2", "R"),
        row("r3", "R"),
        row("r4", "R"),
        row("s1", "S"),
        row("s2", "S"),
        row("s3", "S"),
        row("s4", "S"),
    ]
    variants = [
        variant("r1", "Rv1979c", "p.Arg409Gln"),
        variant("r2", "Rv1979c", "p.Arg409Gln"),
        variant("r3", "Rv1979c", "p.Arg409Gln"),
        variant("s1", "mtrB", "p.Met517Leu"),
    ]
    unsupported = unsupported_variants_by_isolate(variants)
    tradeoff, predictions, selections = run_nested_cv(
        rows,
        ["baseline_pred_R", "baseline_pred_S", "secondary_unsupported"],
        unsupported,
        folds=2,
        repeats=1,
        seed=1,
        min_train_r_carriers=1,
        min_train_fn_carriers=1,
        min_r_fraction=0.05,
        min_odds_ratio=1.0,
        max_p_value=1.0,
    )
    assert len(tradeoff) == 8
    assert len(predictions) == len(rows) * 2
    assert any((row["model"] == NESTED_MODEL for row in tradeoff))
    assert selections


if __name__ == "__main__":
    test_candidate_feature_name_and_columns()
    test_select_train_candidates_excludes_common_background()
    test_run_nested_cv_outputs_predictions_and_selection_rows()
    print("CRyPTIC candidate nested validation tests passed")
