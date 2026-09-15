from __future__ import annotations
import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
from scripts.export_cryptic_fn_mechanism_expansion_audit import (
    COMMON_BACKGROUND_VARIANTS,
    fisher_one_sided_enrichment,
    odds_ratio_haldane,
    variant_key,
)
from scripts.export_cryptic_supervised_upper_bound import (
    TARGET_SPECIFICITIES,
    build_feature_table,
    choose_threshold,
    confusion,
    read_csv,
    ratio,
    score_row,
    selected_features,
    stratified_folds,
    train_bernoulli_nb,
    write_csv,
)

DEFAULT_ISOLATES = Path(
    "reports/publication/external_bdq_phenotype_validation_isolates.csv"
)
DEFAULT_VARIANTS = Path(
    "reports/publication/external_bdq_phenotype_validation_variant_calls.csv"
)
DEFAULT_TRADEOFF_OUTPUT = Path(
    "reports/publication/cryptic_candidate_nested_cv_tradeoff.csv"
)
DEFAULT_PREDICTION_OUTPUT = Path(
    "reports/publication/cryptic_candidate_nested_cv_predictions.csv"
)
DEFAULT_SELECTION_OUTPUT = Path(
    "reports/publication/cryptic_candidate_nested_feature_selection.csv"
)
REFERENCE_MODEL = "reference_no_secondary"
NESTED_MODEL = "nested_candidate_variants"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Evaluate candidate unsupported variants with fold-internal selection."
    )
    parser.add_argument("--isolates", type=Path, default=DEFAULT_ISOLATES)
    parser.add_argument("--variants", type=Path, default=DEFAULT_VARIANTS)
    parser.add_argument("--tradeoff-output", type=Path, default=DEFAULT_TRADEOFF_OUTPUT)
    parser.add_argument(
        "--prediction-output", type=Path, default=DEFAULT_PREDICTION_OUTPUT
    )
    parser.add_argument(
        "--selection-output", type=Path, default=DEFAULT_SELECTION_OUTPUT
    )
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--seed", type=int, default=1729)
    parser.add_argument("--min-train-r-carriers", type=int, default=2)
    parser.add_argument("--min-train-fn-carriers", type=int, default=1)
    parser.add_argument("--min-r-fraction", type=float, default=0.05)
    parser.add_argument("--min-odds-ratio", type=float, default=2.0)
    parser.add_argument("--max-p-value", type=float, default=0.05)
    return parser.parse_args()


def unsupported_variants_by_isolate(
    variants: list[dict[str, str]],
) -> dict[str, set[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    for row in variants:
        if row["call_source"] == "unsupported_variant":
            out[row["unique_id"]].add(variant_key(row))
    return out


def candidate_feature_name(variant: str) -> str:
    return "candidate_variant__" + variant.replace(":", "__").replace(".", "_").replace(
        "-", "_"
    ).replace(">", "_")


def add_candidate_columns(
    rows: list[dict[str, str]],
    candidate_variants: list[str],
    unsupported_by_isolate: dict[str, set[str]],
) -> list[str]:
    feature_names = []
    for candidate in candidate_variants:
        name = candidate_feature_name(candidate)
        feature_names.append(name)
        for row in rows:
            row[name] = (
                "1"
                if candidate in unsupported_by_isolate.get(row["unique_id"], set())
                else "0"
            )
    return feature_names


def select_train_candidates(
    rows: list[dict[str, str]],
    train_indices: list[int],
    unsupported_by_isolate: dict[str, set[str]],
    min_train_r_carriers: int,
    min_train_fn_carriers: int,
    min_r_fraction: float,
    min_odds_ratio: float,
    max_p_value: float,
) -> list[dict[str, str]]:
    train_set = set(train_indices)
    r_indices = {
        idx for idx in train_indices if rows[idx]["bdq_binary_phenotype"] == "R"
    }
    s_indices = train_set - r_indices
    fn_indices = {idx for idx in r_indices if rows[idx]["baseline_prediction"] == "S"}
    variant_to_indices: dict[str, set[int]] = defaultdict(set)
    id_to_index = {rows[idx]["unique_id"]: idx for idx in train_indices}
    for unique_id, variants in unsupported_by_isolate.items():
        idx = id_to_index.get(unique_id)
        if idx is None:
            continue
        for variant in variants:
            variant_to_indices[variant].add(idx)
    selected = []
    for candidate, carrier_indices in variant_to_indices.items():
        if candidate in COMMON_BACKGROUND_VARIANTS:
            continue
        r_carriers = carrier_indices & r_indices
        s_carriers = carrier_indices & s_indices
        fn_carriers = carrier_indices & fn_indices
        a = len(r_carriers)
        b = len(s_carriers)
        c = len(r_indices) - a
        d = len(s_indices) - b
        r_fraction = a / (a + b) if a + b else 0.0
        odds = odds_ratio_haldane(a, b, c, d)
        p_value = fisher_one_sided_enrichment(a, b, c, d)
        if a < min_train_r_carriers:
            continue
        if len(fn_carriers) < min_train_fn_carriers:
            continue
        if r_fraction < min_r_fraction:
            continue
        if odds < min_odds_ratio:
            continue
        if p_value > max_p_value:
            continue
        selected.append(
            {
                "candidate": candidate,
                "train_carriers": str(len(carrier_indices)),
                "train_r_carriers": str(a),
                "train_s_carriers": str(b),
                "train_fn_carriers": str(len(fn_carriers)),
                "train_r_fraction": f"{r_fraction:.4f}",
                "train_odds_ratio_haldane": f"{odds:.6f}",
                "train_fisher_one_sided_p": f"{p_value:.6g}",
            }
        )
    selected.sort(
        key=lambda row: (
            -int(row["train_fn_carriers"]),
            float(row["train_fisher_one_sided_p"]),
            -float(row["train_odds_ratio_haldane"]),
            row["candidate"],
        )
    )
    return selected


def summarize_counts(
    model: str, target: float, counts: Counter[str], folds: int, repeats: int
) -> dict[str, str]:
    total = counts["tp"] + counts["tn"] + counts["fp"] + counts["fn"]
    return {
        "model": model,
        "target_specificity": f"{target:.2f}",
        "folds": str(folds),
        "repeats": str(repeats),
        "total_predictions": str(total),
        "tp": str(counts["tp"]),
        "tn": str(counts["tn"]),
        "fp": str(counts["fp"]),
        "fn": str(counts["fn"]),
        "recall": ratio(counts["tp"], counts["tp"] + counts["fn"]),
        "specificity": ratio(counts["tn"], counts["tn"] + counts["fp"]),
        "precision": ratio(counts["tp"], counts["tp"] + counts["fp"]),
        "accuracy": ratio(counts["tp"] + counts["tn"], total),
    }


def run_nested_cv(
    rows: list[dict[str, str]],
    base_feature_names: list[str],
    unsupported_by_isolate: dict[str, set[str]],
    folds: int,
    repeats: int,
    seed: int,
    min_train_r_carriers: int,
    min_train_fn_carriers: int,
    min_r_fraction: float,
    min_odds_ratio: float,
    max_p_value: float,
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    fold_sets = stratified_folds(rows, folds, repeats, seed)
    reference_features = selected_features(
        base_feature_names, exclude_groups={"secondary"}
    )
    aggregate = {
        REFERENCE_MODEL: {target: Counter() for target in TARGET_SPECIFICITIES},
        NESTED_MODEL: {target: Counter() for target in TARGET_SPECIFICITIES},
    }
    prediction_rows: list[dict[str, str]] = []
    selection_rows: list[dict[str, str]] = []
    for repeat_idx, split in enumerate(fold_sets, start=1):
        for fold_idx, test_indices in enumerate(split, start=1):
            test_set = set(test_indices)
            train_indices = [idx for idx in range(len(rows)) if idx not in test_set]
            selected_candidates = select_train_candidates(
                rows,
                train_indices,
                unsupported_by_isolate,
                min_train_r_carriers,
                min_train_fn_carriers,
                min_r_fraction,
                min_odds_ratio,
                max_p_value,
            )
            candidate_variants = [row["candidate"] for row in selected_candidates]
            candidate_features = add_candidate_columns(
                rows, candidate_variants, unsupported_by_isolate
            )
            for rank, row in enumerate(selected_candidates, start=1):
                selection_rows.append(
                    {
                        "repeat": str(repeat_idx),
                        "fold": str(fold_idx),
                        "rank": str(rank),
                        **row,
                    }
                )
            model_features = {
                REFERENCE_MODEL: reference_features,
                NESTED_MODEL: [*reference_features, *candidate_features],
            }
            for model_name, features in model_features.items():
                model = train_bernoulli_nb(rows, train_indices, features)
                scores = {
                    idx: score_row(rows[idx], model, features)
                    for idx in train_indices + test_indices
                }
                thresholds = {
                    target: choose_threshold(rows, train_indices, scores, target)
                    for target in TARGET_SPECIFICITIES
                }
                for target, threshold in thresholds.items():
                    aggregate[model_name][target].update(
                        confusion(rows, test_indices, scores, threshold)
                    )
                for idx in test_indices:
                    pred_row = {
                        "model": model_name,
                        "repeat": str(repeat_idx),
                        "fold": str(fold_idx),
                        "unique_id": rows[idx]["unique_id"],
                        "bdq_binary_phenotype": rows[idx]["bdq_binary_phenotype"],
                        "score": f"{scores[idx]:.6f}",
                        "selected_candidate_count": str(len(candidate_features)),
                        "selected_candidates": ";".join(candidate_variants),
                    }
                    for target, threshold in thresholds.items():
                        pred_row[f"pred_at_spec_{target:.2f}"] = (
                            "R" if scores[idx] >= threshold else "S"
                        )
                    prediction_rows.append(pred_row)
    tradeoff_rows = []
    for model_name in (REFERENCE_MODEL, NESTED_MODEL):
        for target in TARGET_SPECIFICITIES:
            tradeoff_rows.append(
                summarize_counts(
                    model_name, target, aggregate[model_name][target], folds, repeats
                )
            )
    return (tradeoff_rows, prediction_rows, selection_rows)


def main() -> int:
    args = parse_args()
    isolates = read_csv(args.isolates)
    variants = read_csv(args.variants)
    rows, base_feature_names = build_feature_table(isolates, variants)
    unsupported = unsupported_variants_by_isolate(variants)
    tradeoff_rows, prediction_rows, selection_rows = run_nested_cv(
        rows,
        base_feature_names,
        unsupported,
        args.folds,
        args.repeats,
        args.seed,
        args.min_train_r_carriers,
        args.min_train_fn_carriers,
        args.min_r_fraction,
        args.min_odds_ratio,
        args.max_p_value,
    )
    write_csv(args.tradeoff_output, tradeoff_rows)
    write_csv(args.prediction_output, prediction_rows)
    write_csv(args.selection_output, selection_rows)
    print(f"tradeoff={args.tradeoff_output}")
    print(f"predictions={args.prediction_output}")
    print(f"prediction_rows={len(prediction_rows)}")
    print(f"selection={args.selection_output}")
    print(f"selection_rows={len(selection_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
