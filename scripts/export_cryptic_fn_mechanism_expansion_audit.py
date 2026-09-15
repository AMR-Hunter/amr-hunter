from __future__ import annotations
import argparse
import csv
import math
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path

DEFAULT_ISOLATES = Path(
    "reports/publication/external_bdq_phenotype_validation_isolates.csv"
)
DEFAULT_VARIANTS = Path(
    "reports/publication/external_bdq_phenotype_validation_variant_calls.csv"
)
DEFAULT_FN_MIC_OUTPUT = Path("reports/publication/cryptic_fn_mic_expansion_audit.csv")
DEFAULT_FN_MIC_SUMMARY_OUTPUT = Path(
    "reports/publication/cryptic_fn_mic_expansion_summary.csv"
)
DEFAULT_VARIANT_OUTPUT = Path(
    "reports/publication/cryptic_fn_recurrent_unsupported_variants.csv"
)
DEFAULT_ENRICHMENT_OUTPUT = Path(
    "reports/publication/cryptic_unsupported_variant_enrichment.csv"
)
DEFAULT_COMBO_OUTPUT = Path(
    "reports/publication/cryptic_gene_combination_enrichment.csv"
)
DEFAULT_PRIORITY_OUTPUT = Path(
    "reports/publication/cryptic_recall_candidate_priority.csv"
)
CORE_GENES = {"Rv0678", "atpE", "pepQ"}
SECONDARY_GENES = {"glpK", "Rv1979c", "mmpL5", "mmpS5"}
CONTEXT_GENES = {"mtrA", "mtrB"}
COMMON_BACKGROUND_VARIANTS = {
    "mtrB:p.Met517Leu",
    "mtrB:p.Pro18Ser",
    "glpK:c.1134C>T",
    "glpK:c.510C>T",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Find recurrent unsupported variants and combinations in CRyPTIC BDQ false negatives."
    )
    parser.add_argument("--isolates", type=Path, default=DEFAULT_ISOLATES)
    parser.add_argument("--variants", type=Path, default=DEFAULT_VARIANTS)
    parser.add_argument("--fn-mic-output", type=Path, default=DEFAULT_FN_MIC_OUTPUT)
    parser.add_argument(
        "--fn-mic-summary-output", type=Path, default=DEFAULT_FN_MIC_SUMMARY_OUTPUT
    )
    parser.add_argument("--variant-output", type=Path, default=DEFAULT_VARIANT_OUTPUT)
    parser.add_argument(
        "--enrichment-output", type=Path, default=DEFAULT_ENRICHMENT_OUTPUT
    )
    parser.add_argument("--combo-output", type=Path, default=DEFAULT_COMBO_OUTPUT)
    parser.add_argument("--priority-output", type=Path, default=DEFAULT_PRIORITY_OUTPUT)
    parser.add_argument("--min-fn-carriers", type=int, default=2)
    parser.add_argument("--min-carriers", type=int, default=2)
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


def parse_mic(value: object) -> tuple[float, str]:
    text = str(value or "").strip()
    if not text:
        return (math.nan, "missing")
    qualifier = ""
    if text.startswith((">", "<")):
        qualifier = text[0]
        text = text[1:]
    if text.startswith("="):
        qualifier = "="
        text = text[1:]
    try:
        return (float(text), qualifier)
    except ValueError:
        return (math.nan, "unparsed")


def mic_bin(value: object) -> str:
    mic, qualifier = parse_mic(value)
    if math.isnan(mic):
        return "missing_or_unparsed"
    if qualifier == ">":
        if mic >= 2:
            return ">2"
        if mic >= 1:
            return ">1_to_2"
        return ">0.25_to_1"
    if mic <= 0.25:
        return "<=0.25"
    if mic <= 0.5:
        return "0.5"
    if mic <= 1:
        return "1.0"
    if mic <= 2:
        return "2.0"
    return ">2"


def variant_key(row: dict[str, str]) -> str:
    return f"{row['gene']}:{row['mutation']}"


def is_false_negative(isolate: dict[str, str]) -> bool:
    return (
        isolate["bdq_binary_phenotype"] == "R"
        and isolate["amr_hunter_isolate_prediction"] == "S"
    )


def is_true_positive(isolate: dict[str, str]) -> bool:
    return (
        isolate["bdq_binary_phenotype"] == "R"
        and isolate["amr_hunter_isolate_prediction"] == "R"
    )


def unsupported_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return [row for row in rows if row["call_source"] == "unsupported_variant"]


def log_comb(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def logsumexp(values: list[float]) -> float:
    finite = [value for value in values if value != float("-inf")]
    if not finite:
        return float("-inf")
    max_value = max(finite)
    return max_value + math.log(sum((math.exp(value - max_value) for value in finite)))


def fisher_one_sided_enrichment(a: int, b: int, c: int, d: int) -> float:
    total = a + b + c + d
    total_r = a + c
    carriers = a + b
    upper = min(total_r, carriers)
    logs = [
        log_comb(total_r, observed_r)
        + log_comb(total - total_r, carriers - observed_r)
        - log_comb(total, carriers)
        for observed_r in range(a, upper + 1)
    ]
    return min(1.0, math.exp(logsumexp(logs)))


def odds_ratio_haldane(a: int, b: int, c: int, d: int) -> float:
    return (a + 0.5) * (d + 0.5) / ((b + 0.5) * (c + 0.5))


def ratio(num: int, den: int) -> str:
    return f"{num / den:.4f}" if den else ""


def fmt(value: float) -> str:
    return f"{value:.4g}"


def genes_for_group(rows: list[dict[str, str]], genes: set[str]) -> list[str]:
    return sorted({row["gene"] for row in rows if row["gene"] in genes})


def unsupported_variant_set(rows: list[dict[str, str]]) -> set[str]:
    return {variant_key(row) for row in unsupported_rows(rows)}


def unsupported_gene_set(rows: list[dict[str, str]]) -> set[str]:
    return {row["gene"] for row in unsupported_rows(rows)}


def build_indices(
    isolates: list[dict[str, str]], variants: list[dict[str, str]]
) -> tuple[dict[str, dict[str, str]], dict[str, list[dict[str, str]]]]:
    isolates_by_id = {row["unique_id"]: row for row in isolates}
    variants_by_id: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in variants:
        variants_by_id[row["unique_id"]].append(row)
    return (isolates_by_id, variants_by_id)


def build_fn_mic_rows(
    isolates: list[dict[str, str]], variants_by_id: dict[str, list[dict[str, str]]]
) -> list[dict[str, str]]:
    rows = []
    for isolate in isolates:
        if not is_false_negative(isolate):
            continue
        variants = variants_by_id.get(isolate["unique_id"], [])
        unsupported = unsupported_rows(variants)
        core = genes_for_group(unsupported, CORE_GENES)
        secondary = genes_for_group(unsupported, SECONDARY_GENES)
        context = genes_for_group(unsupported, CONTEXT_GENES)
        labels = sorted(unsupported_variant_set(variants))
        rows.append(
            {
                "unique_id": isolate["unique_id"],
                "bdq_mic": isolate["bdq_mic"],
                "mic_bin": mic_bin(isolate["bdq_mic"]),
                "unsupported_variant_count": str(len(labels)),
                "unsupported_gene_count": str(len(unsupported_gene_set(variants))),
                "core_unsupported_genes": ";".join(core),
                "secondary_unsupported_genes": ";".join(secondary),
                "context_unsupported_genes": ";".join(context),
                "has_mtrB_Met517Leu": str("mtrB:p.Met517Leu" in labels).lower(),
                "has_mtrB_Pro18Ser": str("mtrB:p.Pro18Ser" in labels).lower(),
                "unsupported_variants": ";".join(labels),
            }
        )
    return sorted(rows, key=lambda row: (row["mic_bin"], row["unique_id"]))


def build_fn_mic_summary_rows(
    fn_mic_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    stats: dict[str, Counter[str]] = defaultdict(Counter)
    for row in fn_mic_rows:
        bin_name = row["mic_bin"]
        stats[bin_name]["fn_isolates"] += 1
        if row["core_unsupported_genes"]:
            stats[bin_name]["with_core_unsupported"] += 1
        if row["secondary_unsupported_genes"]:
            stats[bin_name]["with_secondary_unsupported"] += 1
        if row["context_unsupported_genes"]:
            stats[bin_name]["with_context_unsupported"] += 1
        if row["has_mtrB_Met517Leu"] == "true":
            stats[bin_name]["with_mtrB_Met517Leu"] += 1
        if row["has_mtrB_Pro18Ser"] == "true":
            stats[bin_name]["with_mtrB_Pro18Ser"] += 1
    order = {
        "<=0.25": 0,
        "0.5": 1,
        "1.0": 2,
        ">0.25_to_1": 3,
        "2.0": 4,
        ">1_to_2": 5,
        ">2": 6,
        "missing_or_unparsed": 7,
    }
    rows = []
    for bin_name, counts in sorted(
        stats.items(), key=lambda item: order.get(item[0], 99)
    ):
        total = counts["fn_isolates"]
        rows.append(
            {
                "mic_bin": bin_name,
                "fn_isolates": str(total),
                "with_core_unsupported": str(counts["with_core_unsupported"]),
                "with_secondary_unsupported": str(counts["with_secondary_unsupported"]),
                "with_context_unsupported": str(counts["with_context_unsupported"]),
                "with_mtrB_Met517Leu": str(counts["with_mtrB_Met517Leu"]),
                "with_mtrB_Pro18Ser": str(counts["with_mtrB_Pro18Ser"]),
                "fraction_of_fn": ratio(total, len(fn_mic_rows)),
            }
        )
    return rows


def carrier_sets_by_variant(
    variants_by_id: dict[str, list[dict[str, str]]],
) -> dict[str, set[str]]:
    carriers: dict[str, set[str]] = defaultdict(set)
    for unique_id, rows in variants_by_id.items():
        for key in unsupported_variant_set(rows):
            carriers[key].add(unique_id)
    return carriers


def variant_enrichment_rows(
    isolates: list[dict[str, str]],
    variants_by_id: dict[str, list[dict[str, str]]],
    min_carriers: int,
) -> list[dict[str, str]]:
    isolate_ids = {row["unique_id"] for row in isolates}
    r_ids = {row["unique_id"] for row in isolates if row["bdq_binary_phenotype"] == "R"}
    s_ids = {row["unique_id"] for row in isolates if row["bdq_binary_phenotype"] == "S"}
    fn_ids = {row["unique_id"] for row in isolates if is_false_negative(row)}
    tp_ids = {row["unique_id"] for row in isolates if is_true_positive(row)}
    total_r = len(r_ids)
    total_s = len(s_ids)
    rows = []
    for key, raw_carriers in carrier_sets_by_variant(variants_by_id).items():
        carriers = raw_carriers & isolate_ids
        if len(carriers) < min_carriers:
            continue
        r_carriers = carriers & r_ids
        s_carriers = carriers & s_ids
        a = len(r_carriers)
        b = len(s_carriers)
        c = total_r - a
        d = total_s - b
        gene, mutation = key.split(":", 1)
        rows.append(
            {
                "variant": key,
                "gene": gene,
                "mutation": mutation,
                "carrier_total": str(len(carriers)),
                "r_carriers": str(a),
                "s_carriers": str(b),
                "fn_carriers": str(len(carriers & fn_ids)),
                "tp_carriers": str(len(carriers & tp_ids)),
                "r_fraction_among_carriers": ratio(a, a + b),
                "r_recall_gap_fraction": ratio(len(carriers & fn_ids), total_r),
                "odds_ratio_haldane": fmt(odds_ratio_haldane(a, b, c, d)),
                "fisher_one_sided_p": fmt(fisher_one_sided_enrichment(a, b, c, d)),
                "is_common_background": str(key in COMMON_BACKGROUND_VARIANTS).lower(),
            }
        )
    rows.sort(
        key=lambda row: (
            -int(row["fn_carriers"]),
            float(row["fisher_one_sided_p"]),
            -float(row["odds_ratio_haldane"]),
            row["variant"],
        )
    )
    return rows


def recurrent_fn_variant_rows(
    enrichment_rows: list[dict[str, str]], min_fn_carriers: int
) -> list[dict[str, str]]:
    rows = [
        row for row in enrichment_rows if int(row["fn_carriers"]) >= min_fn_carriers
    ]
    return sorted(
        rows,
        key=lambda row: (
            row["is_common_background"] == "true",
            -int(row["fn_carriers"]),
            float(row["fisher_one_sided_p"]),
        ),
    )


def combo_memberships(rows: list[dict[str, str]]) -> set[str]:
    unsupported_genes = unsupported_gene_set(rows)
    unsupported_variants = unsupported_variant_set(rows)
    memberships = set()
    for gene in sorted(unsupported_genes):
        memberships.add(f"gene:{gene}")
    for left, right in combinations(sorted(unsupported_genes), 2):
        memberships.add(f"gene_pair:{left}+{right}")
    if "mtrB:p.Met517Leu" in unsupported_variants:
        memberships.add("context:mtrB_Met517Leu")
    if "mtrB:p.Pro18Ser" in unsupported_variants:
        memberships.add("context:mtrB_Pro18Ser")
    if {"mtrB:p.Met517Leu", "mtrB:p.Pro18Ser"} <= unsupported_variants:
        memberships.add("context:mtrB_double_hit")
    if (
        "mtrB:p.Met517Leu" in unsupported_variants
        and unsupported_genes & SECONDARY_GENES
    ):
        memberships.add("context:mtrB_Met517Leu_plus_secondary_gene")
    if "mtrB:p.Met517Leu" in unsupported_variants and unsupported_genes & CORE_GENES:
        memberships.add("context:mtrB_Met517Leu_plus_core_gene")
    return memberships


def combination_enrichment_rows(
    isolates: list[dict[str, str]],
    variants_by_id: dict[str, list[dict[str, str]]],
    min_carriers: int,
) -> list[dict[str, str]]:
    r_ids = {row["unique_id"] for row in isolates if row["bdq_binary_phenotype"] == "R"}
    s_ids = {row["unique_id"] for row in isolates if row["bdq_binary_phenotype"] == "S"}
    fn_ids = {row["unique_id"] for row in isolates if is_false_negative(row)}
    total_r = len(r_ids)
    total_s = len(s_ids)
    carriers_by_combo: dict[str, set[str]] = defaultdict(set)
    for isolate in isolates:
        unique_id = isolate["unique_id"]
        for combo in combo_memberships(variants_by_id.get(unique_id, [])):
            carriers_by_combo[combo].add(unique_id)
    rows = []
    for combo, carriers in carriers_by_combo.items():
        if len(carriers) < min_carriers:
            continue
        r_carriers = carriers & r_ids
        s_carriers = carriers & s_ids
        a = len(r_carriers)
        b = len(s_carriers)
        c = total_r - a
        d = total_s - b
        rows.append(
            {
                "combination": combo,
                "carrier_total": str(len(carriers)),
                "r_carriers": str(a),
                "s_carriers": str(b),
                "fn_carriers": str(len(carriers & fn_ids)),
                "r_fraction_among_carriers": ratio(a, a + b),
                "r_recall_gap_fraction": ratio(len(carriers & fn_ids), total_r),
                "odds_ratio_haldane": fmt(odds_ratio_haldane(a, b, c, d)),
                "fisher_one_sided_p": fmt(fisher_one_sided_enrichment(a, b, c, d)),
            }
        )
    rows.sort(
        key=lambda row: (
            -int(row["fn_carriers"]),
            float(row["fisher_one_sided_p"]),
            -float(row["odds_ratio_haldane"]),
            row["combination"],
        )
    )
    return rows


def priority_tier(
    fn_carriers: int, r_fraction: float, odds_ratio: float, is_common: bool
) -> str:
    if is_common and r_fraction < 0.05:
        return "background_context_only"
    if fn_carriers >= 2 and r_fraction >= 0.2 and (odds_ratio >= 5):
        return "high_priority_candidate"
    if fn_carriers >= 2 and r_fraction >= 0.05 and (odds_ratio >= 2):
        return "medium_priority_candidate"
    if fn_carriers >= 2:
        return "recurrent_low_specificity_context"
    return "low_priority"


def priority_rows(
    variant_rows: list[dict[str, str]], combo_rows: list[dict[str, str]]
) -> list[dict[str, str]]:
    rows = []
    for source, source_rows, name_field in (
        ("unsupported_variant", variant_rows, "variant"),
        ("gene_or_context_combination", combo_rows, "combination"),
    ):
        for row in source_rows:
            fn_carriers = int(row["fn_carriers"])
            if fn_carriers < 2:
                continue
            r_fraction = float(row["r_fraction_among_carriers"] or 0)
            odds = float(row["odds_ratio_haldane"])
            is_common = row.get("is_common_background", "false") == "true"
            tier = priority_tier(fn_carriers, r_fraction, odds, is_common)
            if tier == "low_priority":
                continue
            rows.append(
                {
                    "candidate": row[name_field],
                    "candidate_type": source,
                    "priority_tier": tier,
                    "fn_carriers": row["fn_carriers"],
                    "r_carriers": row["r_carriers"],
                    "s_carriers": row["s_carriers"],
                    "carrier_total": row["carrier_total"],
                    "r_fraction_among_carriers": row["r_fraction_among_carriers"],
                    "odds_ratio_haldane": row["odds_ratio_haldane"],
                    "fisher_one_sided_p": row["fisher_one_sided_p"],
                    "recommended_next_step": recommended_next_step(tier, source),
                }
            )
    order = {
        "high_priority_candidate": 0,
        "medium_priority_candidate": 1,
        "recurrent_low_specificity_context": 2,
        "background_context_only": 3,
    }
    rows.sort(
        key=lambda row: (
            order[row["priority_tier"]],
            -int(row["fn_carriers"]),
            float(row["fisher_one_sided_p"]),
            row["candidate"],
        )
    )
    return rows


def recommended_next_step(tier: str, source: str) -> str:
    if tier == "high_priority_candidate":
        return "literature_review_and_targeted_modeling_then_nested_validation"
    if tier == "medium_priority_candidate":
        return "inspect_carrier_isolates_and_add_only_as_nested_cv_feature"
    if tier == "recurrent_low_specificity_context":
        return "keep_as_context_or_indeterminate_signal_not_resistant_rule"
    if source == "unsupported_variant":
        return "treat_as_background_until_external_evidence"
    return "use_only_for_context_stratification"


def main() -> int:
    args = parse_args()
    isolates = read_csv(args.isolates)
    variants = read_csv(args.variants)
    _isolates_by_id, variants_by_id = build_indices(isolates, variants)
    fn_mic_rows = build_fn_mic_rows(isolates, variants_by_id)
    fn_mic_summary_rows = build_fn_mic_summary_rows(fn_mic_rows)
    enrichment_rows = variant_enrichment_rows(
        isolates, variants_by_id, args.min_carriers
    )
    recurrent_rows = recurrent_fn_variant_rows(enrichment_rows, args.min_fn_carriers)
    combo_rows = combination_enrichment_rows(
        isolates, variants_by_id, args.min_carriers
    )
    candidate_rows = priority_rows(recurrent_rows, combo_rows)
    write_csv(args.fn_mic_output, fn_mic_rows)
    write_csv(args.fn_mic_summary_output, fn_mic_summary_rows)
    write_csv(args.enrichment_output, enrichment_rows)
    write_csv(args.variant_output, recurrent_rows)
    write_csv(args.combo_output, combo_rows)
    write_csv(args.priority_output, candidate_rows)
    print(f"fn_mic={args.fn_mic_output}")
    print(f"fn_rows={len(fn_mic_rows)}")
    print(f"fn_mic_summary={args.fn_mic_summary_output}")
    print(f"fn_mic_summary_rows={len(fn_mic_summary_rows)}")
    print(f"unsupported_variant_enrichment={args.enrichment_output}")
    print(f"variant_rows={len(enrichment_rows)}")
    print(f"recurrent_fn_variants={args.variant_output}")
    print(f"recurrent_rows={len(recurrent_rows)}")
    print(f"combination_enrichment={args.combo_output}")
    print(f"combo_rows={len(combo_rows)}")
    print(f"candidate_priority={args.priority_output}")
    print(f"candidate_rows={len(candidate_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
