from __future__ import annotations
import argparse
import csv
import math
import random
from collections import Counter, defaultdict
from pathlib import Path

DEFAULT_ISOLATES = Path(
    "reports/publication/external_bdq_phenotype_validation_isolates.csv"
)
DEFAULT_VARIANTS = Path(
    "reports/publication/external_bdq_phenotype_validation_variant_calls.csv"
)
DEFAULT_FEATURE_OUTPUT = Path("reports/publication/cryptic_isolate_feature_table.csv")
DEFAULT_TRADEOFF_OUTPUT = Path("reports/publication/cryptic_supervised_cv_tradeoff.csv")
DEFAULT_PREDICTION_OUTPUT = Path(
    "reports/publication/cryptic_supervised_cv_predictions.csv"
)
DEFAULT_WEIGHT_OUTPUT = Path(
    "reports/publication/cryptic_supervised_feature_weights.csv"
)
CORE_GENES = {"Rv0678", "atpE", "pepQ"}
SECONDARY_GENES = {"glpK", "Rv1979c", "mmpL5", "mmpS5", "lpqB"}
CONTEXT_GENES = {"mtrA", "mtrB"}
ALL_GENES = sorted(CORE_GENES | SECONDARY_GENES | CONTEXT_GENES)
assert len(ALL_GENES) == 10, f"expected the 10 target genes, got {ALL_GENES}"
TARGET_SPECIFICITIES = (0.99, 0.98, 0.95, 0.9)
MODEL_SPECS = {
    "all_features": {"include": None, "exclude": set()},
    "baseline_only": {"include": {"baseline"}, "exclude": set()},
    "no_mtrb": {"include": None, "exclude": {"mtrb"}},
    "no_secondary": {"include": None, "exclude": {"secondary"}},
    "no_core": {"include": None, "exclude": {"core"}},
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build genotype-only isolate features and run supervised CRyPTIC CV."
    )
    parser.add_argument("--isolates", type=Path, default=DEFAULT_ISOLATES)
    parser.add_argument("--variants", type=Path, default=DEFAULT_VARIANTS)
    parser.add_argument("--feature-output", type=Path, default=DEFAULT_FEATURE_OUTPUT)
    parser.add_argument("--tradeoff-output", type=Path, default=DEFAULT_TRADEOFF_OUTPUT)
    parser.add_argument(
        "--prediction-output", type=Path, default=DEFAULT_PREDICTION_OUTPUT
    )
    parser.add_argument("--weight-output", type=Path, default=DEFAULT_WEIGHT_OUTPUT)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--seed", type=int, default=1729)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(
    path: Path, rows: list[dict[str, str]], fieldnames: list[str] | None = None
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def add_feature(features: dict[str, int], name: str, value: bool | int) -> None:
    features[name] = int(bool(value))


def count_where(rows: list[dict[str, str]], **conditions: str) -> int:
    count = 0
    for row in rows:
        if all((row.get(key) == value for key, value in conditions.items())):
            count += 1
    return count


def gene_rows(rows: list[dict[str, str]], genes: set[str]) -> list[dict[str, str]]:
    return [row for row in rows if row["gene"] in genes]


def make_features(
    rows: list[dict[str, str]], baseline_prediction: str
) -> dict[str, int]:
    features: dict[str, int] = {}
    variant_count = len(rows)
    unsupported = [row for row in rows if row["call_source"] == "unsupported_variant"]
    r_calls = [row for row in rows if row["amr_hunter_prediction"] == "R"]
    s_calls = [row for row in rows if row["amr_hunter_prediction"] == "S"]
    core_rows = gene_rows(rows, CORE_GENES)
    secondary_rows = gene_rows(rows, SECONDARY_GENES)
    context_rows = gene_rows(rows, CONTEXT_GENES)
    unsupported_core = [row for row in unsupported if row["gene"] in CORE_GENES]
    unsupported_secondary = [
        row for row in unsupported if row["gene"] in SECONDARY_GENES
    ]
    unsupported_context = [row for row in unsupported if row["gene"] in CONTEXT_GENES]
    add_feature(features, "baseline_pred_R", baseline_prediction == "R")
    add_feature(features, "baseline_pred_S", baseline_prediction == "S")
    add_feature(features, "baseline_pred_UNKNOWN", baseline_prediction == "UNKNOWN")
    add_feature(features, "variant_count_ge_2", variant_count >= 2)
    add_feature(features, "variant_count_ge_4", variant_count >= 4)
    add_feature(features, "variant_count_ge_6", variant_count >= 6)
    add_feature(features, "unsupported_count_ge_1", len(unsupported) >= 1)
    add_feature(features, "unsupported_count_ge_2", len(unsupported) >= 2)
    add_feature(features, "unsupported_count_ge_3", len(unsupported) >= 3)
    add_feature(features, "r_call_count_ge_1", len(r_calls) >= 1)
    add_feature(features, "r_call_count_ge_2", len(r_calls) >= 2)
    add_feature(features, "s_call_count_ge_1", len(s_calls) >= 1)
    add_feature(
        features,
        "direct_s_count_ge_1",
        count_where(rows, call_source="direct_record_match", amr_hunter_prediction="S")
        >= 1,
    )
    add_feature(
        features,
        "snapshot_r_count_ge_1",
        count_where(rows, call_source="snapshot_virtual", amr_hunter_prediction="R")
        >= 1,
    )
    add_feature(
        features,
        "consequence_r_count_ge_1",
        count_where(rows, call_source="consequence_inferred", amr_hunter_prediction="R")
        >= 1,
    )
    add_feature(
        features,
        "direct_r_count_ge_1",
        count_where(rows, call_source="direct_record_match", amr_hunter_prediction="R")
        >= 1,
    )
    add_feature(features, "core_present", bool(core_rows))
    add_feature(features, "secondary_present", bool(secondary_rows))
    add_feature(features, "context_present", bool(context_rows))
    add_feature(features, "core_unsupported", bool(unsupported_core))
    add_feature(features, "secondary_unsupported", bool(unsupported_secondary))
    add_feature(features, "context_unsupported", bool(unsupported_context))
    add_feature(
        features,
        "core_r_call",
        any((row["amr_hunter_prediction"] == "R" for row in core_rows)),
    )
    add_feature(
        features,
        "secondary_r_call",
        any((row["amr_hunter_prediction"] == "R" for row in secondary_rows)),
    )
    add_feature(
        features,
        "context_r_call",
        any((row["amr_hunter_prediction"] == "R" for row in context_rows)),
    )
    keys = {(row["gene"], row["mutation"]) for row in unsupported}
    add_feature(features, "mtrB_Met517Leu", ("mtrB", "p.Met517Leu") in keys)
    add_feature(features, "mtrB_Pro18Ser", ("mtrB", "p.Pro18Ser") in keys)
    add_feature(
        features,
        "mtrB_double_hit",
        ("mtrB", "p.Met517Leu") in keys and ("mtrB", "p.Pro18Ser") in keys,
    )
    add_feature(
        features,
        "mtrB_Met517Leu_plus_core",
        ("mtrB", "p.Met517Leu") in keys and bool(unsupported_core),
    )
    add_feature(
        features,
        "mtrB_Met517Leu_plus_secondary",
        ("mtrB", "p.Met517Leu") in keys and bool(unsupported_secondary),
    )
    for gene in ALL_GENES:
        rows_gene = [row for row in rows if row["gene"] == gene]
        add_feature(features, f"gene_{gene}_present", bool(rows_gene))
        add_feature(
            features,
            f"gene_{gene}_unsupported",
            any((row["call_source"] == "unsupported_variant" for row in rows_gene)),
        )
        add_feature(
            features,
            f"gene_{gene}_r_call",
            any((row["amr_hunter_prediction"] == "R" for row in rows_gene)),
        )
        add_feature(
            features,
            f"gene_{gene}_frameshift",
            any((row["variant_class"] == "frameshift" for row in rows_gene)),
        )
    return features


def feature_group(name: str) -> str:
    if name.startswith("baseline_"):
        return "baseline"
    if name.startswith(
        (
            "variant_",
            "unsupported_count",
            "r_call_count",
            "s_call_count",
            "direct_",
            "snapshot_",
            "consequence_",
        )
    ):
        return "burden"
    if name.startswith("mtrB_"):
        return "mtrb"
    if name.startswith("gene_"):
        for gene in CORE_GENES:
            if name.startswith(f"gene_{gene}_"):
                return "core"
        for gene in SECONDARY_GENES:
            if name.startswith(f"gene_{gene}_"):
                return "secondary"
        return "context"
    if name.startswith("core_"):
        return "core"
    if name.startswith("secondary_"):
        return "secondary"
    if name.startswith("context_"):
        return "context"
    return "other"


def build_feature_table(
    isolates: list[dict[str, str]], variants: list[dict[str, str]]
) -> tuple[list[dict[str, str]], list[str]]:
    variants_by_isolate: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in variants:
        variants_by_isolate[row["unique_id"]].append(row)
    feature_dicts: list[dict[str, int]] = []
    for isolate in isolates:
        feature_dicts.append(
            make_features(
                variants_by_isolate.get(isolate["unique_id"], []),
                isolate["amr_hunter_isolate_prediction"],
            )
        )
    feature_names = sorted({key for features in feature_dicts for key in features})
    rows_out: list[dict[str, str]] = []
    for isolate, features in zip(isolates, feature_dicts):
        out = {
            "unique_id": isolate["unique_id"],
            "ena_sample": isolate["ena_sample"],
            "bdq_binary_phenotype": isolate["bdq_binary_phenotype"],
            "label": "1" if isolate["bdq_binary_phenotype"] == "R" else "0",
            "bdq_mic": isolate["bdq_mic"],
            "baseline_prediction": isolate["amr_hunter_isolate_prediction"],
        }
        for name in feature_names:
            out[name] = str(features.get(name, 0))
        rows_out.append(out)
    return (rows_out, feature_names)


def stratified_folds(
    rows: list[dict[str, str]], folds: int, repeats: int, seed: int
) -> list[list[list[int]]]:
    positives = [index for index, row in enumerate(rows) if row["label"] == "1"]
    negatives = [index for index, row in enumerate(rows) if row["label"] == "0"]
    all_repeats: list[list[list[int]]] = []
    for repeat in range(repeats):
        rng = random.Random(seed + repeat)
        pos = positives[:]
        neg = negatives[:]
        rng.shuffle(pos)
        rng.shuffle(neg)
        split = [[] for _ in range(folds)]
        for idx, row_index in enumerate(pos):
            split[idx % folds].append(row_index)
        for idx, row_index in enumerate(neg):
            split[idx % folds].append(row_index)
        for fold in split:
            rng.shuffle(fold)
        all_repeats.append(split)
    return all_repeats


def selected_features(
    feature_names: list[str],
    include_groups: set[str] | None = None,
    exclude_groups: set[str] | None = None,
) -> list[str]:
    excluded = exclude_groups or set()
    selected = []
    for name in feature_names:
        group = feature_group(name)
        if include_groups is not None and group not in include_groups:
            continue
        if group in excluded:
            continue
        selected.append(name)
    return selected


def train_bernoulli_nb(
    rows: list[dict[str, str]], indices: list[int], features: list[str]
) -> dict[str, float]:
    pos_indices = [idx for idx in indices if rows[idx]["label"] == "1"]
    neg_indices = [idx for idx in indices if rows[idx]["label"] == "0"]
    pos_n = len(pos_indices)
    neg_n = len(neg_indices)
    model: dict[str, float] = {"__prior__": math.log((pos_n + 1) / (neg_n + 1))}
    for feature in features:
        pos_present = sum((1 for idx in pos_indices if rows[idx][feature] == "1"))
        neg_present = sum((1 for idx in neg_indices if rows[idx][feature] == "1"))
        p_pos = (pos_present + 1) / (pos_n + 2)
        p_neg = (neg_present + 1) / (neg_n + 2)
        present_weight = math.log(p_pos / p_neg)
        absent_weight = math.log((1 - p_pos) / (1 - p_neg))
        model[feature] = present_weight
        model[f"!{feature}"] = absent_weight
    return model


def score_row(
    row: dict[str, str], model: dict[str, float], features: list[str]
) -> float:
    score = model["__prior__"]
    for feature in features:
        score += model[feature] if row[feature] == "1" else model[f"!{feature}"]
    return score


def confusion(
    rows: list[dict[str, str]],
    indices: list[int],
    scores: dict[int, float],
    threshold: float,
) -> Counter[str]:
    out: Counter[str] = Counter()
    for idx in indices:
        truth = rows[idx]["label"] == "1"
        pred = scores[idx] >= threshold
        if truth and pred:
            out["tp"] += 1
        elif truth and (not pred):
            out["fn"] += 1
        elif not truth and pred:
            out["fp"] += 1
        else:
            out["tn"] += 1
    return out


def ratio(num: int, den: int) -> str:
    return f"{num / den:.4f}" if den else ""


def specificity(counter: Counter[str]) -> float:
    den = counter["tn"] + counter["fp"]
    return counter["tn"] / den if den else 0.0


def recall(counter: Counter[str]) -> float:
    den = counter["tp"] + counter["fn"]
    return counter["tp"] / den if den else 0.0


def choose_threshold(
    rows: list[dict[str, str]],
    train_indices: list[int],
    scores: dict[int, float],
    target_specificity: float,
) -> float:
    total_pos = sum((1 for idx in train_indices if rows[idx]["label"] == "1"))
    total_neg = len(train_indices) - total_pos
    best_threshold = float("inf")
    best_recall = -1.0
    best_specificity = -1.0

    def consider(threshold: float, tp: int, fp: int) -> None:
        nonlocal best_threshold, best_recall, best_specificity
        fn = total_pos - tp
        tn = total_neg - fp
        counts = Counter({"tp": tp, "fp": fp, "fn": fn, "tn": tn})
        spec = specificity(counts)
        rec = recall(counts)
        if spec + 1e-12 < target_specificity:
            return
        if rec > best_recall or (rec == best_recall and spec > best_specificity):
            best_threshold = threshold
            best_recall = rec
            best_specificity = spec

    tp = 0
    fp = 0
    consider(float("inf"), tp, fp)
    sorted_indices = sorted(train_indices, key=lambda idx: scores[idx], reverse=True)
    position = 0
    while position < len(sorted_indices):
        threshold = scores[sorted_indices[position]]
        while (
            position < len(sorted_indices)
            and scores[sorted_indices[position]] == threshold
        ):
            idx = sorted_indices[position]
            if rows[idx]["label"] == "1":
                tp += 1
            else:
                fp += 1
            position += 1
        consider(threshold, tp, fp)
    consider(float("-inf"), tp, fp)
    return best_threshold


def run_cv(
    rows: list[dict[str, str]],
    feature_names: list[str],
    folds: int,
    repeats: int,
    seed: int,
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    fold_sets = stratified_folds(rows, folds, repeats, seed)
    tradeoff_rows: list[dict[str, str]] = []
    prediction_rows: list[dict[str, str]] = []
    weight_accumulator: dict[str, Counter[str]] = defaultdict(Counter)
    for model_name, spec in MODEL_SPECS.items():
        features = selected_features(
            feature_names,
            include_groups=spec["include"],
            exclude_groups=spec["exclude"],
        )
        aggregate: dict[float, Counter[str]] = {
            target: Counter() for target in TARGET_SPECIFICITIES
        }
        thresholds_by_target: dict[float, list[float]] = {
            target: [] for target in TARGET_SPECIFICITIES
        }
        for repeat_idx, split in enumerate(fold_sets, start=1):
            for fold_idx, test_indices in enumerate(split, start=1):
                test_set = set(test_indices)
                train_indices = [idx for idx in range(len(rows)) if idx not in test_set]
                model = train_bernoulli_nb(rows, train_indices, features)
                for feature in features:
                    weight_accumulator[model_name][feature] += model[feature]
                    weight_accumulator[model_name][f"__n__:{feature}"] += 1
                scores = {
                    idx: score_row(rows[idx], model, features)
                    for idx in train_indices + test_indices
                }
                for target in TARGET_SPECIFICITIES:
                    threshold = choose_threshold(rows, train_indices, scores, target)
                    thresholds_by_target[target].append(threshold)
                    counts = confusion(rows, test_indices, scores, threshold)
                    aggregate[target].update(counts)
                for idx in test_indices:
                    pred_row = {
                        "model": model_name,
                        "repeat": str(repeat_idx),
                        "fold": str(fold_idx),
                        "unique_id": rows[idx]["unique_id"],
                        "bdq_binary_phenotype": rows[idx]["bdq_binary_phenotype"],
                        "score": f"{scores[idx]:.6f}",
                    }
                    for target in TARGET_SPECIFICITIES:
                        threshold = thresholds_by_target[target][-1]
                        pred_row[f"pred_at_spec_{target:.2f}"] = (
                            "R" if scores[idx] >= threshold else "S"
                        )
                    prediction_rows.append(pred_row)
        for target in TARGET_SPECIFICITIES:
            counts = aggregate[target]
            total = counts["tp"] + counts["tn"] + counts["fp"] + counts["fn"]
            tradeoff_rows.append(
                {
                    "model": model_name,
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
                    "mean_threshold": f"{sum(thresholds_by_target[target]) / len(thresholds_by_target[target]):.6f}",
                }
            )
    weight_rows: list[dict[str, str]] = []
    for model_name, weights in weight_accumulator.items():
        for feature in feature_names:
            n = weights[f"__n__:{feature}"]
            if n == 0:
                continue
            weight_rows.append(
                {
                    "model": model_name,
                    "feature": feature,
                    "feature_group": feature_group(feature),
                    "mean_present_log_odds_weight": f"{weights[feature] / n:.6f}",
                }
            )
    weight_rows.sort(
        key=lambda row: (row["model"], -abs(float(row["mean_present_log_odds_weight"])))
    )
    return (tradeoff_rows, prediction_rows, weight_rows)


def main() -> int:
    args = parse_args()
    if args.folds < 2:
        raise ValueError("--folds must be at least 2")
    if args.repeats < 1:
        raise ValueError("--repeats must be at least 1")
    isolates = read_csv(args.isolates)
    variants = read_csv(args.variants)
    feature_rows, feature_names = build_feature_table(isolates, variants)
    feature_fieldnames = [
        "unique_id",
        "ena_sample",
        "bdq_binary_phenotype",
        "label",
        "bdq_mic",
        "baseline_prediction",
        *feature_names,
    ]
    write_csv(args.feature_output, feature_rows, feature_fieldnames)
    tradeoff_rows, prediction_rows, weight_rows = run_cv(
        feature_rows, feature_names, args.folds, args.repeats, args.seed
    )
    write_csv(args.tradeoff_output, tradeoff_rows)
    write_csv(args.prediction_output, prediction_rows)
    write_csv(args.weight_output, weight_rows)
    print(f"features={args.feature_output}")
    print(f"feature_rows={len(feature_rows)}")
    print(f"feature_count={len(feature_names)}")
    print(f"tradeoff={args.tradeoff_output}")
    print(f"prediction_rows={len(prediction_rows)}")
    print(f"weights={args.weight_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
