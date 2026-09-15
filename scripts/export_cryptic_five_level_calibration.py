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
DEFAULT_ISOLATE_OUTPUT = Path("reports/publication/cryptic_five_level_isolates.csv")
DEFAULT_LABEL_OUTPUT = Path("reports/publication/cryptic_five_level_label_summary.csv")
DEFAULT_METRIC_OUTPUT = Path("reports/publication/cryptic_five_level_metrics.csv")
CORE_GENES = {"Rv0678", "atpE", "pepQ"}
SECONDARY_GENES = {"glpK", "Rv1979c", "mmpL5", "mmpS5"}
CONTEXT_GENES = {"mtrA", "mtrB"}
LABEL_ORDER = ["R", "Likely R", "Indeterminate", "Likely S", "S"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export five-level CRyPTIC BDQ isolate calibration labels."
    )
    parser.add_argument("--isolates", type=Path, default=DEFAULT_ISOLATES)
    parser.add_argument("--variants", type=Path, default=DEFAULT_VARIANTS)
    parser.add_argument("--isolate-output", type=Path, default=DEFAULT_ISOLATE_OUTPUT)
    parser.add_argument("--label-output", type=Path, default=DEFAULT_LABEL_OUTPUT)
    parser.add_argument("--metric-output", type=Path, default=DEFAULT_METRIC_OUTPUT)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def variant_key(row: dict[str, str]) -> tuple[str, str]:
    return (row["gene"], row["mutation"])


def variant_label(
    row: dict[str, str], carrier_stats: Counter[str] | None = None
) -> str:
    base = f"{row['gene']}:{row['mutation']}"
    if carrier_stats is None:
        return base
    total = carrier_stats["R"] + carrier_stats["S"]
    return f"{base}({carrier_stats['R']}/{total})"


def carrier_fraction(stats: Counter[str]) -> float:
    total = stats["R"] + stats["S"]
    return stats["R"] / total if total else 0.0


def split_resistant_evidence(
    rows: list[dict[str, str]], carrier_stats: dict[tuple[str, str], Counter[str]]
) -> tuple[list[str], list[str], list[str]]:
    strong: list[str] = []
    likely: list[str] = []
    weak: list[str] = []
    for row in rows:
        if row["amr_hunter_prediction"] != "R":
            continue
        stats = carrier_stats[variant_key(row)]
        r_count = stats["R"]
        fraction = carrier_fraction(stats)
        label = variant_label(row, stats)
        if r_count >= 1 and fraction >= 0.25:
            strong.append(label)
        elif r_count >= 1 and (
            fraction >= 0.1
            or (
                row["call_source"] == "direct_record_match"
                and row["gene"] in CORE_GENES
            )
        ):
            likely.append(label)
        else:
            weak.append(label)
    return (strong, likely, weak)


def split_unsupported_evidence(
    rows: list[dict[str, str]], carrier_stats: dict[tuple[str, str], Counter[str]]
) -> tuple[list[str], list[str], list[str]]:
    likely_r: list[str] = []
    indeterminate: list[str] = []
    common_context: list[str] = []
    for row in rows:
        if row["call_source"] != "unsupported_variant":
            continue
        stats = carrier_stats[variant_key(row)]
        r_count = stats["R"]
        s_count = stats["S"]
        fraction = carrier_fraction(stats)
        label = variant_label(row, stats)
        if (
            row["gene"] not in CONTEXT_GENES
            and r_count >= 1
            and (s_count <= 3 or fraction >= 0.25)
        ):
            likely_r.append(label)
        elif row["gene"] in CORE_GENES:
            indeterminate.append(label)
        else:
            common_context.append(label)
    return (likely_r, indeterminate, common_context)


def assign_five_level_label(
    isolate: dict[str, str],
    rows: list[dict[str, str]],
    carrier_stats: dict[tuple[str, str], Counter[str]],
) -> dict[str, str]:
    baseline = isolate["amr_hunter_isolate_prediction"]
    strong_r, likely_r_from_r_call, weak_r = split_resistant_evidence(
        rows, carrier_stats
    )
    likely_r_from_unsupported, core_unsupported, common_context = (
        split_unsupported_evidence(rows, carrier_stats)
    )
    unsupported_count = sum(
        (1 for row in rows if row["call_source"] == "unsupported_variant")
    )
    direct_sensitive = [
        variant_label(row)
        for row in rows
        if row["call_source"] == "direct_record_match"
        and row["amr_hunter_prediction"] == "S"
    ]
    if baseline == "UNKNOWN":
        label = "Indeterminate"
        reason = "baseline_unknown"
        evidence = ""
    elif strong_r:
        label = "R"
        reason = "carrier_enriched_resistant_signal"
        evidence = ";".join(strong_r)
    elif likely_r_from_r_call or likely_r_from_unsupported:
        label = "Likely R"
        reason = "moderate_or_rare_resistant_signal"
        evidence = ";".join((likely_r_from_r_call + likely_r_from_unsupported)[:5])
    elif weak_r:
        label = "Indeterminate"
        reason = "weak_or_conflicting_resistant_signal"
        evidence = ";".join(weak_r[:5])
    elif core_unsupported:
        label = "Indeterminate"
        reason = "unsupported_core_gene_signal"
        evidence = ";".join(core_unsupported[:5])
    elif unsupported_count:
        label = "Likely S"
        reason = "common_or_context_unsupported_without_resistant_signal"
        evidence = ";".join(common_context[:5])
    else:
        label = "S"
        reason = "sensitive_markers_without_unsupported_signal"
        evidence = ";".join(direct_sensitive[:5])
    return {
        "dataset": isolate["dataset"],
        "unique_id": isolate["unique_id"],
        "ena_sample": isolate["ena_sample"],
        "bdq_binary_phenotype": isolate["bdq_binary_phenotype"],
        "bdq_mic": isolate["bdq_mic"],
        "baseline_prediction": baseline,
        "five_level_label": label,
        "label_reason": reason,
        "label_evidence": evidence,
        "variant_count": str(len(rows)),
        "unsupported_count": str(unsupported_count),
        "strong_r_evidence": ";".join(strong_r),
        "likely_r_evidence": ";".join(
            (likely_r_from_r_call + likely_r_from_unsupported)[:10]
        ),
        "weak_r_evidence": ";".join(weak_r[:10]),
        "core_unsupported_evidence": ";".join(core_unsupported[:10]),
        "common_context_evidence": ";".join(common_context[:10]),
        "baseline_call_source_summary": isolate["call_source_summary"],
    }


def build_carrier_stats(
    variants: list[dict[str, str]],
) -> dict[tuple[str, str], Counter[str]]:
    stats: dict[tuple[str, str], Counter[str]] = defaultdict(Counter)
    for row in variants:
        stats[variant_key(row)][row["bdq_binary_phenotype"]] += 1
    return stats


def build_label_summary(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    total = len(rows)
    phenotype_r = sum((1 for row in rows if row["bdq_binary_phenotype"] == "R"))
    phenotype_s = sum((1 for row in rows if row["bdq_binary_phenotype"] == "S"))
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    for row in rows:
        counts[row["five_level_label"]][row["bdq_binary_phenotype"]] += 1
    out: list[dict[str, str]] = []
    for label in LABEL_ORDER:
        r_count = counts[label]["R"]
        s_count = counts[label]["S"]
        label_total = r_count + s_count
        out.append(
            {
                "label": label,
                "total_isolates": str(label_total),
                "phenotype_r": str(r_count),
                "phenotype_s": str(s_count),
                "fraction_of_all": f"{label_total / total:.4f}" if total else "",
                "fraction_of_resistant": (
                    f"{r_count / phenotype_r:.4f}" if phenotype_r else ""
                ),
                "fraction_of_sensitive": (
                    f"{s_count / phenotype_s:.4f}" if phenotype_s else ""
                ),
                "r_fraction_within_label": (
                    f"{r_count / label_total:.4f}" if label_total else ""
                ),
            }
        )
    return out


def build_metric_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    for row in rows:
        counts[row["five_level_label"]][row["bdq_binary_phenotype"]] += 1
    phenotype_r = sum((counts[label]["R"] for label in LABEL_ORDER))
    phenotype_s = sum((counts[label]["S"] for label in LABEL_ORDER))

    def count(labels: set[str], phenotype: str) -> int:
        return sum((counts[label][phenotype] for label in labels))

    def ratio(numerator: int, denominator: int) -> str:
        return f"{numerator / denominator:.4f}" if denominator else ""

    metrics = [
        (
            "hard_r_precision",
            count({"R"}, "R"),
            count({"R"}, "R") + count({"R"}, "S"),
            "Phenotype-R fraction among hard R labels.",
        ),
        (
            "likely_r_precision",
            count({"Likely R"}, "R"),
            count({"Likely R"}, "R") + count({"Likely R"}, "S"),
            "Phenotype-R fraction among Likely R labels.",
        ),
        (
            "r_or_likely_r_recall",
            count({"R", "Likely R"}, "R"),
            phenotype_r,
            "Resistant isolates captured as R or Likely R.",
        ),
        (
            "r_not_called_sensitive",
            count({"R", "Likely R", "Indeterminate"}, "R"),
            phenotype_r,
            "Resistant isolates not assigned hard/likely sensitive labels.",
        ),
        (
            "s_not_called_resistant_or_likely",
            count({"Indeterminate", "Likely S", "S"}, "S"),
            phenotype_s,
            "Sensitive isolates not assigned R or Likely R.",
        ),
        (
            "indeterminate_rate",
            count({"Indeterminate"}, "R") + count({"Indeterminate"}, "S"),
            phenotype_r + phenotype_s,
            "Fraction of cohort assigned Indeterminate.",
        ),
        (
            "hard_s_risk",
            count({"S"}, "R"),
            count({"S"}, "R") + count({"S"}, "S"),
            "Phenotype-R fraction among hard S labels.",
        ),
        (
            "likely_s_risk",
            count({"Likely S"}, "R"),
            count({"Likely S"}, "R") + count({"Likely S"}, "S"),
            "Phenotype-R fraction among Likely S labels.",
        ),
    ]
    return [
        {
            "metric": name,
            "numerator": str(numerator),
            "denominator": str(denominator),
            "value": ratio(numerator, denominator),
            "notes": notes,
        }
        for name, numerator, denominator, notes in metrics
    ]


def main() -> int:
    args = parse_args()
    isolates = read_csv(args.isolates)
    variants = read_csv(args.variants)
    variants_by_isolate: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in variants:
        variants_by_isolate[row["unique_id"]].append(row)
    carrier_stats = build_carrier_stats(variants)
    isolate_rows = [
        assign_five_level_label(
            isolate, variants_by_isolate.get(isolate["unique_id"], []), carrier_stats
        )
        for isolate in isolates
    ]
    write_csv(args.isolate_output, isolate_rows)
    write_csv(args.label_output, build_label_summary(isolate_rows))
    write_csv(args.metric_output, build_metric_rows(isolate_rows))
    print(f"isolates={args.isolate_output}")
    print(f"isolate_rows={len(isolate_rows)}")
    print(f"label_summary={args.label_output}")
    print(f"metrics={args.metric_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
