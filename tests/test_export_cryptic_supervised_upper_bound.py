from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
from scripts.export_cryptic_supervised_upper_bound import (
    build_feature_table,
    choose_threshold,
    confusion,
    feature_group,
    recall,
    score_row,
    selected_features,
    specificity,
    stratified_folds,
    train_bernoulli_nb,
)


def isolate(unique_id: str, phenotype: str, prediction: str) -> dict[str, str]:
    return {
        "unique_id": unique_id,
        "ena_sample": unique_id,
        "bdq_binary_phenotype": phenotype,
        "bdq_mic": "0.5" if phenotype == "R" else "0.06",
        "amr_hunter_isolate_prediction": prediction,
    }


def variant(
    unique_id: str,
    gene: str,
    mutation: str,
    prediction: str = "UNKNOWN",
    source: str = "unsupported_variant",
    variant_class: str = "missense",
) -> dict[str, str]:
    return {
        "unique_id": unique_id,
        "gene": gene,
        "mutation": mutation,
        "amr_hunter_prediction": prediction,
        "call_source": source,
        "variant_class": variant_class,
    }


def test_build_feature_table_extracts_genotype_only_signals():
    isolates = [isolate("iso-r", "R", "R"), isolate("iso-s", "S", "S")]
    variants = [
        variant("iso-r", "mtrB", "p.Met517Leu"),
        variant("iso-r", "mtrB", "p.Pro18Ser"),
        variant("iso-r", "Rv0678", "p.Arg90Cys", "R", "snapshot_virtual"),
        variant("iso-r", "Rv0678", "p.Leu10fs", variant_class="frameshift"),
        variant("iso-s", "glpK", "p.Val460Ala"),
    ]
    rows, features = build_feature_table(isolates, variants)
    by_id = {row["unique_id"]: row for row in rows}
    assert "mtrB_double_hit" in features
    assert by_id["iso-r"]["baseline_pred_R"] == "1"
    assert by_id["iso-r"]["mtrB_Met517Leu"] == "1"
    assert by_id["iso-r"]["mtrB_double_hit"] == "1"
    assert by_id["iso-r"]["core_unsupported"] == "1"
    assert by_id["iso-r"]["gene_Rv0678_frameshift"] == "1"
    assert by_id["iso-r"]["snapshot_r_count_ge_1"] == "1"
    assert by_id["iso-s"]["secondary_unsupported"] == "1"
    assert by_id["iso-s"]["mtrB_double_hit"] == "0"


def test_selected_features_supports_baseline_only_and_ablation_specs():
    feature_names = [
        "baseline_pred_R",
        "mtrB_Met517Leu",
        "secondary_unsupported",
        "gene_Rv0678_unsupported",
        "variant_count_ge_4",
    ]
    baseline = selected_features(feature_names, include_groups={"baseline"})
    no_mtrb = selected_features(feature_names, exclude_groups={"mtrb"})
    assert baseline == ["baseline_pred_R"]
    assert "mtrB_Met517Leu" not in no_mtrb
    assert "secondary_unsupported" in no_mtrb
    assert feature_group("gene_Rv0678_unsupported") == "core"


def test_cross_validation_helpers_choose_specificity_constrained_threshold():
    rows = [
        {"unique_id": "r1", "label": "1", "signal": "1", "bdq_binary_phenotype": "R"},
        {"unique_id": "r2", "label": "1", "signal": "1", "bdq_binary_phenotype": "R"},
        {"unique_id": "s1", "label": "0", "signal": "0", "bdq_binary_phenotype": "S"},
        {"unique_id": "s2", "label": "0", "signal": "0", "bdq_binary_phenotype": "S"},
    ]
    indices = list(range(len(rows)))
    model = train_bernoulli_nb(rows, indices, ["signal"])
    scores = {idx: score_row(rows[idx], model, ["signal"]) for idx in indices}
    threshold = choose_threshold(rows, indices, scores, target_specificity=1.0)
    counts = confusion(rows, indices, scores, threshold)
    assert specificity(counts) == 1.0
    assert recall(counts) == 1.0


def test_stratified_folds_preserve_all_indices_once_per_repeat():
    rows = [
        {"label": "1"},
        {"label": "1"},
        {"label": "0"},
        {"label": "0"},
        {"label": "0"},
        {"label": "0"},
    ]
    folds = stratified_folds(rows, folds=2, repeats=2, seed=7)
    assert len(folds) == 2
    for split in folds:
        seen = sorted((index for fold in split for index in fold))
        assert seen == list(range(len(rows)))


if __name__ == "__main__":
    test_build_feature_table_extracts_genotype_only_signals()
    test_selected_features_supports_baseline_only_and_ablation_specs()
    test_cross_validation_helpers_choose_specificity_constrained_threshold()
    test_stratified_folds_preserve_all_indices_once_per_repeat()
    print("CRyPTIC supervised upper-bound tests passed")
