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
DEFAULT_FIVE_LEVEL = Path("reports/publication/cryptic_five_level_isolates.csv")
DEFAULT_GROUP_OUTPUT = Path("reports/publication/cryptic_mtrb_enrichment_summary.csv")
DEFAULT_COMBO_OUTPUT = Path("reports/publication/cryptic_mtrb_combo_summary.csv")
DEFAULT_DETAIL_OUTPUT = Path("reports/publication/cryptic_mtrb_context_isolates.csv")
CORE_GENES = {"Rv0678", "atpE", "pepQ"}
SECONDARY_GENES = {"glpK", "Rv1979c", "mmpL5", "mmpS5"}
MTRB_MET517 = ("mtrB", "p.Met517Leu")
MTRB_PRO18 = ("mtrB", "p.Pro18Ser")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build mtrB/mtrAB enrichment summaries from CRyPTIC BDQ validation outputs."
    )
    parser.add_argument("--isolates", type=Path, default=DEFAULT_ISOLATES)
    parser.add_argument("--variants", type=Path, default=DEFAULT_VARIANTS)
    parser.add_argument("--five-level", type=Path, default=DEFAULT_FIVE_LEVEL)
    parser.add_argument("--group-output", type=Path, default=DEFAULT_GROUP_OUTPUT)
    parser.add_argument("--combo-output", type=Path, default=DEFAULT_COMBO_OUTPUT)
    parser.add_argument("--detail-output", type=Path, default=DEFAULT_DETAIL_OUTPUT)
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


def variant_key(row: dict[str, str]) -> tuple[str, str]:
    return (row["gene"], row["mutation"])


def unsupported_keys(rows: list[dict[str, str]]) -> set[tuple[str, str]]:
    return {
        variant_key(row) for row in rows if row["call_source"] == "unsupported_variant"
    }


def mtrb_signature(rows: list[dict[str, str]]) -> str:
    labels = sorted(
        (
            f"{row['mutation']}[{row['amr_hunter_prediction']}/{row['call_source']}]"
            for row in rows
            if row["gene"] == "mtrB"
        )
    )
    return " | ".join(labels) if labels else "(no_mtrB_variant)"


def has_unsupported_gene(rows: list[dict[str, str]], genes: set[str]) -> bool:
    return any(
        (
            row["call_source"] == "unsupported_variant" and row["gene"] in genes
            for row in rows
        )
    )


def group_memberships(rows: list[dict[str, str]]) -> set[str]:
    keys = unsupported_keys(rows)
    has_met517 = MTRB_MET517 in keys
    has_pro18 = MTRB_PRO18 in keys
    has_core = has_unsupported_gene(rows, CORE_GENES)
    has_secondary = has_unsupported_gene(rows, SECONDARY_GENES)
    groups = {"all"}
    if has_met517 or has_pro18:
        groups.add("mtrB_context_any")
    else:
        groups.add("mtrB_context_absent")
    if has_met517:
        groups.add("mtrB_Met517Leu")
    if has_pro18:
        groups.add("mtrB_Pro18Ser")
    if has_met517 and has_pro18:
        groups.add("mtrB_double_hit")
    if has_met517 and (not has_pro18):
        groups.add("mtrB_Met517Leu_only")
    if has_pro18 and (not has_met517):
        groups.add("mtrB_Pro18Ser_only")
    if has_met517 and has_core:
        groups.add("mtrB_Met517Leu_plus_core_unsupported")
    if has_met517 and has_secondary:
        groups.add("mtrB_Met517Leu_plus_secondary_unsupported")
    if has_met517 and has_pro18 and has_secondary:
        groups.add("mtrB_double_hit_plus_secondary_unsupported")
    return groups


def join_counter(counter: Counter[str]) -> str:
    return ";".join((f"{key}:{value}" for key, value in sorted(counter.items())))


def summarize_group(
    name: str,
    isolate_rows: list[dict[str, str]],
    five_level_by_id: dict[str, dict[str, str]],
) -> dict[str, str]:
    total = len(isolate_rows)
    phenotype = Counter((row["bdq_binary_phenotype"] for row in isolate_rows))
    baseline = Counter((row["amr_hunter_isolate_prediction"] for row in isolate_rows))
    five = Counter(
        (
            five_level_by_id.get(row["unique_id"], {}).get(
                "five_level_label", "missing"
            )
            for row in isolate_rows
        )
    )
    mic_all = Counter((mic_bin(row["bdq_mic"]) for row in isolate_rows))
    mic_r = Counter(
        (
            mic_bin(row["bdq_mic"])
            for row in isolate_rows
            if row["bdq_binary_phenotype"] == "R"
        )
    )
    r_fraction = phenotype["R"] / total if total else 0.0
    return {
        "group": name,
        "total_isolates": str(total),
        "phenotype_r": str(phenotype["R"]),
        "phenotype_s": str(phenotype["S"]),
        "r_fraction": f"{r_fraction:.4f}" if total else "",
        "baseline_pred_R": str(baseline["R"]),
        "baseline_pred_S": str(baseline["S"]),
        "baseline_pred_UNKNOWN": str(baseline["UNKNOWN"]),
        "five_level_R": str(five["R"]),
        "five_level_Likely_R": str(five["Likely R"]),
        "five_level_Indeterminate": str(five["Indeterminate"]),
        "five_level_Likely_S": str(five["Likely S"]),
        "five_level_S": str(five["S"]),
        "mic_bins_all": join_counter(mic_all),
        "mic_bins_r": join_counter(mic_r),
    }


def build_group_summary(
    isolates: list[dict[str, str]],
    variants_by_isolate: dict[str, list[dict[str, str]]],
    five_level_by_id: dict[str, dict[str, str]],
) -> list[dict[str, str]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for isolate in isolates:
        for group in group_memberships(
            variants_by_isolate.get(isolate["unique_id"], [])
        ):
            grouped[group].append(isolate)
    order = [
        "all",
        "mtrB_context_any",
        "mtrB_context_absent",
        "mtrB_Met517Leu",
        "mtrB_Pro18Ser",
        "mtrB_double_hit",
        "mtrB_Met517Leu_only",
        "mtrB_Pro18Ser_only",
        "mtrB_Met517Leu_plus_core_unsupported",
        "mtrB_Met517Leu_plus_secondary_unsupported",
        "mtrB_double_hit_plus_secondary_unsupported",
    ]
    return [
        summarize_group(name, grouped.get(name, []), five_level_by_id)
        for name in order
        if name in grouped
    ]


def build_combo_summary(
    isolates: list[dict[str, str]],
    variants_by_isolate: dict[str, list[dict[str, str]]],
    five_level_by_id: dict[str, dict[str, str]],
) -> list[dict[str, str]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    secondary_counts: dict[str, Counter[str]] = defaultdict(Counter)
    core_counts: dict[str, Counter[str]] = defaultdict(Counter)
    for isolate in isolates:
        rows = variants_by_isolate.get(isolate["unique_id"], [])
        signature = mtrb_signature(rows)
        grouped[signature].append(isolate)
        for row in rows:
            if row["call_source"] != "unsupported_variant":
                continue
            label = f"{row['gene']}:{row['mutation']}"
            if row["gene"] in SECONDARY_GENES:
                secondary_counts[signature][label] += 1
            if row["gene"] in CORE_GENES:
                core_counts[signature][label] += 1
    combo_rows: list[dict[str, str]] = []
    for signature, rows in grouped.items():
        summary = summarize_group(signature, rows, five_level_by_id)
        summary["mtrB_signature"] = summary.pop("group")
        summary["top_secondary_unsupported"] = ";".join(
            (
                f"{key}:{value}"
                for key, value in secondary_counts[signature].most_common(5)
            )
        )
        summary["top_core_unsupported"] = ";".join(
            (f"{key}:{value}" for key, value in core_counts[signature].most_common(5))
        )
        combo_rows.append(summary)
    combo_rows.sort(
        key=lambda row: (int(row["phenotype_r"]), int(row["total_isolates"])),
        reverse=True,
    )
    return combo_rows


def build_detail_rows(
    isolates: list[dict[str, str]],
    variants_by_isolate: dict[str, list[dict[str, str]]],
    five_level_by_id: dict[str, dict[str, str]],
) -> list[dict[str, str]]:
    rows_out: list[dict[str, str]] = []
    for isolate in isolates:
        rows = variants_by_isolate.get(isolate["unique_id"], [])
        keys = unsupported_keys(rows)
        include = (
            isolate["bdq_binary_phenotype"] == "R"
            or MTRB_PRO18 in keys
            or (
                MTRB_MET517 in keys
                and has_unsupported_gene(rows, CORE_GENES | SECONDARY_GENES)
            )
        )
        if not include:
            continue
        five = five_level_by_id.get(isolate["unique_id"], {})
        rows_out.append(
            {
                "unique_id": isolate["unique_id"],
                "ena_sample": isolate["ena_sample"],
                "bdq_binary_phenotype": isolate["bdq_binary_phenotype"],
                "bdq_mic": isolate["bdq_mic"],
                "mic_bin": mic_bin(isolate["bdq_mic"]),
                "baseline_prediction": isolate["amr_hunter_isolate_prediction"],
                "five_level_label": five.get("five_level_label", ""),
                "mtrB_signature": mtrb_signature(rows),
                "group_memberships": ";".join(sorted(group_memberships(rows))),
                "core_unsupported": ";".join(
                    sorted(
                        (
                            f"{row['gene']}:{row['mutation']}"
                            for row in rows
                            if row["call_source"] == "unsupported_variant"
                            and row["gene"] in CORE_GENES
                        )
                    )
                ),
                "secondary_unsupported": ";".join(
                    sorted(
                        (
                            f"{row['gene']}:{row['mutation']}"
                            for row in rows
                            if row["call_source"] == "unsupported_variant"
                            and row["gene"] in SECONDARY_GENES
                        )
                    )
                ),
                "baseline_call_source_summary": isolate["call_source_summary"],
            }
        )
    return rows_out


def main() -> int:
    args = parse_args()
    isolates = read_csv(args.isolates)
    variants = read_csv(args.variants)
    five_level = read_csv(args.five_level) if args.five_level.exists() else []
    five_level_by_id = {row["unique_id"]: row for row in five_level}
    variants_by_isolate: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in variants:
        variants_by_isolate[row["unique_id"]].append(row)
    group_rows = build_group_summary(isolates, variants_by_isolate, five_level_by_id)
    combo_rows = build_combo_summary(isolates, variants_by_isolate, five_level_by_id)
    detail_rows = build_detail_rows(isolates, variants_by_isolate, five_level_by_id)
    write_csv(args.group_output, group_rows)
    write_csv(args.combo_output, combo_rows)
    write_csv(args.detail_output, detail_rows)
    print(f"group_summary={args.group_output}")
    print(f"group_rows={len(group_rows)}")
    print(f"combo_summary={args.combo_output}")
    print(f"combo_rows={len(combo_rows)}")
    print(f"context_isolates={args.detail_output}")
    print(f"context_rows={len(detail_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
