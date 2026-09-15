from __future__ import annotations
import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

DEFAULT_ISOLATES = Path(
    "reports/publication/external_bdq_phenotype_validation_isolates.csv"
)
DEFAULT_VARIANTS = Path(
    "reports/publication/external_bdq_phenotype_validation_variant_calls.csv"
)
DEFAULT_VARIANT_OUTPUT = Path(
    "reports/publication/cryptic_targeted_gpu_variant_runlist.csv"
)
DEFAULT_COMBO_OUTPUT = Path(
    "reports/publication/cryptic_targeted_gpu_combo_runlist.csv"
)
CORE_GENES = {"Rv0678", "atpE", "pepQ"}
SECONDARY_GENES = {"glpK", "Rv1979c", "mmpL5", "mmpS5"}
DEFERRED_CONTEXT_GENES = {"mtrB", "mtrA"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build targeted GPU runlists from CRyPTIC BDQ validation errors."
    )
    parser.add_argument("--isolates", type=Path, default=DEFAULT_ISOLATES)
    parser.add_argument("--variants", type=Path, default=DEFAULT_VARIANTS)
    parser.add_argument("--variant-output", type=Path, default=DEFAULT_VARIANT_OUTPUT)
    parser.add_argument("--combo-output", type=Path, default=DEFAULT_COMBO_OUTPUT)
    parser.add_argument("--max-variants", type=int, default=80)
    parser.add_argument("--max-combos", type=int, default=80)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def variant_key(row: dict[str, str]) -> tuple[str, str]:
    return (row["gene"], row["mutation"])


def gene_priority(gene: str) -> int:
    if gene in CORE_GENES:
        return 4
    if gene in SECONDARY_GENES:
        return 3
    if gene in DEFERRED_CONTEXT_GENES:
        return 1
    return 2


def source_priority(source: str) -> int:
    if source == "unsupported_variant":
        return 4
    if source == "snapshot_virtual":
        return 3
    if source == "consequence_inferred":
        return 3
    if source == "direct_record_match":
        return 2
    return 1


def variant_class_priority(variant_class: str, region_type: str) -> int:
    text = f"{variant_class} {region_type}".lower()
    if "frameshift" in text or "start_lost" in text or "nonsense" in text:
        return 4
    if "protein_substitution" in text:
        return 3
    if "regulatory" in text:
        return 2
    return 1


def classify_variant(row: dict[str, str], fn_count: int, fp_count: int) -> str:
    gene = row["gene"]
    prediction = row["amr_hunter_prediction"]
    source = row["call_source"]
    if fn_count and source == "unsupported_variant" and (gene in CORE_GENES):
        return "P0_FN_CORE_UNSUPPORTED"
    if (
        fp_count
        and prediction == "R"
        and (source in {"snapshot_virtual", "consequence_inferred"})
    ):
        return "P0_FP_R_CALL_CALIBRATION"
    if fn_count and source == "unsupported_variant" and (gene in SECONDARY_GENES):
        return "P1_FN_SECONDARY_UNSUPPORTED"
    if fn_count and source == "unsupported_variant":
        return "P2_FN_CONTEXT_UNSUPPORTED"
    if fp_count and prediction == "R":
        return "P1_FP_R_CALL_REVIEW"
    return "P3_CONTEXT"


def build_variant_rows(
    isolates: dict[str, dict[str, str]], variants: list[dict[str, str]], max_rows: int
) -> list[dict[str, str]]:
    stats: dict[tuple[str, str], Counter] = defaultdict(Counter)
    exemplar: dict[tuple[str, str], dict[str, str]] = {}
    for row in variants:
        key = variant_key(row)
        isolate = isolates[row["unique_id"]]
        phenotype = isolate["bdq_binary_phenotype"]
        prediction = isolate["amr_hunter_isolate_prediction"]
        stats[key][f"phenotype_{phenotype}"] += 1
        if phenotype == "R" and prediction == "S":
            stats[key]["fn"] += 1
        if phenotype == "S" and prediction == "R":
            stats[key]["fp"] += 1
        if phenotype == "R" and prediction == "R":
            stats[key]["tp"] += 1
        if phenotype == "S" and prediction == "S":
            stats[key]["tn"] += 1
        if row["amr_hunter_prediction"] == "R":
            stats[key]["variant_r_calls"] += 1
        if row["amr_hunter_prediction"] == "S":
            stats[key]["variant_s_calls"] += 1
        if row["amr_hunter_prediction"] == "UNKNOWN":
            stats[key]["variant_unknown_calls"] += 1
        exemplar.setdefault(key, row)
    ranked: list[tuple[tuple[int, int, int, int, int, int], dict[str, str]]] = []
    for key, counts in stats.items():
        row = exemplar[key]
        fn_count = counts["fn"]
        fp_count = counts["fp"]
        if not fn_count and (not fp_count):
            continue
        phenotype_r = counts["phenotype_R"]
        phenotype_s = counts["phenotype_S"]
        total = phenotype_r + phenotype_s
        r_enrichment = phenotype_r / total if total else 0.0
        category = classify_variant(row, fn_count, fp_count)
        p0 = 1 if category.startswith("P0") else 0
        score = (
            p0,
            gene_priority(row["gene"]),
            source_priority(row["call_source"]),
            variant_class_priority(row["variant_class"], row["region_type"]),
            fn_count + fp_count,
            -phenotype_s if fn_count and (not fp_count) else fp_count,
        )
        ranked.append(
            (
                score,
                {
                    "category": category,
                    "gene": row["gene"],
                    "mutation": row["mutation"],
                    "region_type": row["region_type"],
                    "variant_class": row["variant_class"],
                    "consequence": row.get("consequence", ""),
                    "current_variant_prediction": row["amr_hunter_prediction"],
                    "call_source": row["call_source"],
                    "phenotype_r_carriers": str(phenotype_r),
                    "phenotype_s_carriers": str(phenotype_s),
                    "r_carrier_fraction": f"{r_enrichment:.4f}",
                    "fn_isolates": str(fn_count),
                    "fp_isolates": str(fp_count),
                    "tp_isolates": str(counts["tp"]),
                    "tn_isolates": str(counts["tn"]),
                    "variant_r_calls": str(counts["variant_r_calls"]),
                    "variant_s_calls": str(counts["variant_s_calls"]),
                    "variant_unknown_calls": str(counts["variant_unknown_calls"]),
                },
            )
        )
    ranked.sort(key=lambda item: item[0], reverse=True)
    rows = []
    for rank, (_, row) in enumerate(ranked[:max_rows], start=1):
        row = {"rank": str(rank), **row}
        rows.append(row)
    return rows


def combo_signature(rows: list[dict[str, str]]) -> str:
    parts = []
    for row in rows:
        if (
            row["call_source"] != "unsupported_variant"
            and row["amr_hunter_prediction"] == "S"
        ):
            continue
        parts.append(
            ":".join(
                [
                    row["gene"],
                    row["mutation"],
                    row["amr_hunter_prediction"],
                    row["call_source"],
                ]
            )
        )
    if not parts:
        return "(only_current_S_or_common_markers)"
    return " | ".join(sorted(parts))


def build_combo_rows(
    isolates: dict[str, dict[str, str]],
    variants_by_isolate: dict[str, list[dict[str, str]]],
    max_rows: int,
) -> list[dict[str, str]]:
    combos: dict[str, Counter] = defaultdict(Counter)
    examples: dict[str, dict[str, str]] = {}
    genes_by_combo: dict[str, set[str]] = defaultdict(set)
    for unique_id, isolate in isolates.items():
        phenotype = isolate["bdq_binary_phenotype"]
        prediction = isolate["amr_hunter_isolate_prediction"]
        if (phenotype, prediction) not in {("R", "S"), ("S", "R")}:
            continue
        rows = variants_by_isolate.get(unique_id, [])
        signature = combo_signature(rows)
        combos[signature][f"{phenotype}_{prediction}"] += 1
        combos[signature]["isolates"] += 1
        examples.setdefault(signature, isolate)
        for row in rows:
            genes_by_combo[signature].add(row["gene"])
    ranked: list[tuple[tuple[int, int, int], dict[str, str]]] = []
    for signature, counts in combos.items():
        genes = genes_by_combo[signature]
        fn = counts["R_S"]
        fp = counts["S_R"]
        priority = 3 if genes & CORE_GENES else 2 if genes & SECONDARY_GENES else 1
        score = (priority, fn + fp, fn)
        example = examples[signature]
        ranked.append(
            (
                score,
                {
                    "combo_signature": signature,
                    "genes": ";".join(sorted(genes)),
                    "fn_isolates": str(fn),
                    "fp_isolates": str(fp),
                    "error_isolates": str(counts["isolates"]),
                    "example_unique_id": example["unique_id"],
                    "example_bdq_mic": example["bdq_mic"],
                    "example_prediction": example["amr_hunter_isolate_prediction"],
                    "example_call_source_summary": example["call_source_summary"],
                },
            )
        )
    ranked.sort(key=lambda item: item[0], reverse=True)
    rows = []
    for rank, (_, row) in enumerate(ranked[:max_rows], start=1):
        rows.append({"rank": str(rank), **row})
    return rows


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    isolate_rows = read_csv(args.isolates)
    variants = read_csv(args.variants)
    isolates = {row["unique_id"]: row for row in isolate_rows}
    variants_by_isolate: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in variants:
        if row["unique_id"] in isolates:
            variants_by_isolate[row["unique_id"]].append(row)
    variant_rows = build_variant_rows(isolates, variants, args.max_variants)
    combo_rows = build_combo_rows(isolates, variants_by_isolate, args.max_combos)
    write_csv(args.variant_output, variant_rows)
    write_csv(args.combo_output, combo_rows)
    print(f"variant_runlist={args.variant_output}")
    print(f"variant_rows={len(variant_rows)}")
    print(f"combo_runlist={args.combo_output}")
    print(f"combo_rows={len(combo_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
