from __future__ import annotations
import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path

DEFAULT_ISOLATES = Path(
    "reports/publication/external_bdq_phenotype_validation_isolates.csv"
)
DEFAULT_VARIANTS = Path(
    "reports/publication/external_bdq_phenotype_validation_variant_calls.csv"
)
DEFAULT_FN_OUTPUT = Path("reports/publication/cryptic_fn_taxonomy.csv")
DEFAULT_MIC_OUTPUT = Path("reports/publication/cryptic_mic_stratification.csv")
DEFAULT_CATEGORY_OUTPUT = Path("reports/publication/cryptic_error_taxonomy_summary.csv")
CORE_GENES = {"Rv0678", "atpE", "pepQ"}
SECONDARY_GENES = {"glpK", "Rv1979c", "mmpL5", "mmpS5"}
CONTEXT_GENES = {"mtrA", "mtrB"}
COMMON_BACKGROUND_VARIANTS = {
    ("mtrB", "p.Met517Leu"),
    ("mtrB", "p.Pro18Ser"),
    ("glpK", "c.1134C>T"),
    ("glpK", "c.510C>T"),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build CRyPTIC BDQ false-negative taxonomy and MIC summaries."
    )
    parser.add_argument("--isolates", type=Path, default=DEFAULT_ISOLATES)
    parser.add_argument("--variants", type=Path, default=DEFAULT_VARIANTS)
    parser.add_argument("--fn-output", type=Path, default=DEFAULT_FN_OUTPUT)
    parser.add_argument("--mic-output", type=Path, default=DEFAULT_MIC_OUTPUT)
    parser.add_argument("--category-output", type=Path, default=DEFAULT_CATEGORY_OUTPUT)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


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


def variant_label(row: dict[str, str]) -> str:
    return f"{row['gene']}:{row['mutation']}[{row['amr_hunter_prediction']}/{row['call_source']}]"


def unsupported_labels(
    rows: list[dict[str, str]], genes: set[str] | None = None
) -> list[str]:
    labels = []
    for row in rows:
        if row["call_source"] != "unsupported_variant":
            continue
        if genes is not None and row["gene"] not in genes:
            continue
        labels.append(f"{row['gene']}:{row['mutation']}")
    return sorted(labels)


def classify_fn(rows: list[dict[str, str]]) -> str:
    core = unsupported_labels(rows, CORE_GENES)
    secondary = unsupported_labels(rows, SECONDARY_GENES)
    context = unsupported_labels(rows, CONTEXT_GENES)
    unsupported = unsupported_labels(rows)
    has_mtrb_met517 = any(
        (row["gene"] == "mtrB" and row["mutation"] == "p.Met517Leu" for row in rows)
    )
    has_mtrb_pro18 = any(
        (row["gene"] == "mtrB" and row["mutation"] == "p.Pro18Ser" for row in rows)
    )
    has_only_s_calls = rows and all(
        (row["amr_hunter_prediction"] == "S" for row in rows)
    )
    if core:
        return "core_gene_unsupported"
    if has_mtrb_met517 and has_mtrb_pro18:
        return "mtrB_double_hit_context"
    if has_mtrb_met517 and (not secondary):
        return "mtrB_met517_context_only"
    if secondary:
        return "secondary_gene_unsupported"
    if context:
        return "context_gene_unsupported"
    if unsupported:
        return "other_unsupported"
    if has_only_s_calls:
        return "only_current_sensitive_markers"
    return "other"


def build_fn_rows(
    isolates: list[dict[str, str]], variants_by_isolate: dict[str, list[dict[str, str]]]
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for isolate in isolates:
        if (
            isolate["bdq_binary_phenotype"] != "R"
            or isolate["amr_hunter_isolate_prediction"] != "S"
        ):
            continue
        variants = variants_by_isolate.get(isolate["unique_id"], [])
        unsupported = [
            row for row in variants if row["call_source"] == "unsupported_variant"
        ]
        direct_s = [
            row
            for row in variants
            if row["call_source"] == "direct_record_match"
            and row["amr_hunter_prediction"] == "S"
        ]
        rows.append(
            {
                "unique_id": isolate["unique_id"],
                "ena_sample": isolate["ena_sample"],
                "bdq_mic": isolate["bdq_mic"],
                "mic_bin": mic_bin(isolate["bdq_mic"]),
                "taxonomy": classify_fn(variants),
                "variant_count": str(len(variants)),
                "unsupported_count": str(len(unsupported)),
                "direct_sensitive_count": str(len(direct_s)),
                "core_unsupported": ";".join(unsupported_labels(variants, CORE_GENES)),
                "secondary_unsupported": ";".join(
                    unsupported_labels(variants, SECONDARY_GENES)
                ),
                "context_unsupported": ";".join(
                    unsupported_labels(variants, CONTEXT_GENES)
                ),
                "has_mtrB_Met517Leu": str(
                    any(
                        (
                            row["gene"] == "mtrB" and row["mutation"] == "p.Met517Leu"
                            for row in variants
                        )
                    )
                ).lower(),
                "has_mtrB_Pro18Ser": str(
                    any(
                        (
                            row["gene"] == "mtrB" and row["mutation"] == "p.Pro18Ser"
                            for row in variants
                        )
                    )
                ).lower(),
                "common_background_variants": ";".join(
                    sorted(
                        (
                            f"{row['gene']}:{row['mutation']}"
                            for row in variants
                            if (row["gene"], row["mutation"])
                            in COMMON_BACKGROUND_VARIANTS
                        )
                    )
                ),
                "all_interpreted_variants": ";".join(
                    (variant_label(row) for row in variants)
                ),
            }
        )
    return sorted(
        rows, key=lambda row: (row["taxonomy"], row["mic_bin"], row["unique_id"])
    )


def build_mic_rows(isolates: list[dict[str, str]]) -> list[dict[str, str]]:
    stats: dict[tuple[str, str], Counter] = defaultdict(Counter)
    for row in isolates:
        key = (row["bdq_binary_phenotype"], mic_bin(row["bdq_mic"]))
        prediction = row["amr_hunter_isolate_prediction"]
        stats[key]["total"] += 1
        stats[key][f"pred_{prediction}"] += 1
        if row["correct"].lower() == "true":
            stats[key]["correct"] += 1
    rows = []
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
    for (phenotype, bin_name), counts in sorted(
        stats.items(), key=lambda item: (item[0][0], order.get(item[0][1], 99))
    ):
        total = counts["total"]
        if phenotype == "R":
            key_metric = counts["pred_R"] / total if total else 0.0
            metric_name = "sensitivity"
        else:
            key_metric = counts["pred_S"] / total if total else 0.0
            metric_name = "specificity"
        rows.append(
            {
                "phenotype": phenotype,
                "mic_bin": bin_name,
                "total_isolates": str(total),
                "pred_R": str(counts["pred_R"]),
                "pred_S": str(counts["pred_S"]),
                "pred_UNKNOWN": str(counts["pred_UNKNOWN"]),
                "correct": str(counts["correct"]),
                "accuracy": f"{counts['correct'] / total:.4f}" if total else "",
                "primary_metric": metric_name,
                "primary_metric_value": f"{key_metric:.4f}",
            }
        )
    return rows


def build_category_rows(
    fn_rows: list[dict[str, str]], isolates: list[dict[str, str]]
) -> list[dict[str, str]]:
    category_counts = Counter((row["taxonomy"] for row in fn_rows))
    mic_counts: dict[str, Counter] = defaultdict(Counter)
    for row in fn_rows:
        mic_counts[row["taxonomy"]][row["mic_bin"]] += 1
    r_total = sum((1 for row in isolates if row["bdq_binary_phenotype"] == "R"))
    fn_total = len(fn_rows)
    rows = []
    for category, count in category_counts.most_common():
        rows.append(
            {
                "taxonomy": category,
                "fn_isolates": str(count),
                "fraction_of_fn": f"{count / fn_total:.4f}" if fn_total else "",
                "fraction_of_all_r": f"{count / r_total:.4f}" if r_total else "",
                "mic_bin_counts": ";".join(
                    (
                        f"{name}:{value}"
                        for name, value in sorted(mic_counts[category].items())
                    )
                ),
            }
        )
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
    isolates = read_csv(args.isolates)
    variants = read_csv(args.variants)
    variants_by_isolate: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in variants:
        variants_by_isolate[row["unique_id"]].append(row)
    fn_rows = build_fn_rows(isolates, variants_by_isolate)
    mic_rows = build_mic_rows(isolates)
    category_rows = build_category_rows(fn_rows, isolates)
    write_csv(args.fn_output, fn_rows)
    write_csv(args.mic_output, mic_rows)
    write_csv(args.category_output, category_rows)
    print(f"fn_taxonomy={args.fn_output}")
    print(f"fn_rows={len(fn_rows)}")
    print(f"mic_stratification={args.mic_output}")
    print(f"mic_rows={len(mic_rows)}")
    print(f"category_summary={args.category_output}")
    print(f"category_rows={len(category_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
