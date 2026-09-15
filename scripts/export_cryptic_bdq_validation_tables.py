from __future__ import annotations
import argparse
from collections import Counter, defaultdict
from pathlib import Path
import re
import sys
from typing import Iterable, Mapping

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
from core.consequence import resolve_consequence_level_decision
from scripts.cryptic_bdq_validation_common import (
    BDQ_BINARY_PHENOTYPE_FIELD,
    BDQ_MIC_FIELD,
    BDQ_PHENOTYPE_QUALITY_FIELD,
    DATASET_NAME,
    DEFAULT_CRYPTIC_METADATA,
    DEFAULT_VALIDATION_ISOLATES_OUTPUT,
    DEFAULT_VALIDATION_SUMMARY_OUTPUT,
    DEFAULT_VALIDATION_VARIANTS_OUTPUT,
    ENA_SAMPLE_FIELD,
    REQUIRED_METADATA_FIELDS,
    UNIQUE_ID_FIELD,
    VCF_FIELD,
    ValidationMetrics,
    is_present,
    load_bdq_gene_logic,
    normalize_bdq_phenotype,
    normalize_quality_label,
    passes_bdq_quality,
    read_csv_rows,
    read_delimited_rows,
    require_fields,
    write_csv_rows,
)

DEFAULT_CONFIG = Path("config/config.yaml")
DEFAULT_WHO_UNIFIED = Path("reports/publication/who_bdq_per_row_benchmark_unified.csv")
DEFAULT_SNAPSHOT = Path("reports/amr_results_snapshot.csv")
AA_ONE_TO_THREE = {
    "A": "Ala",
    "C": "Cys",
    "D": "Asp",
    "E": "Glu",
    "F": "Phe",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "K": "Lys",
    "L": "Leu",
    "M": "Met",
    "N": "Asn",
    "P": "Pro",
    "Q": "Gln",
    "R": "Arg",
    "S": "Ser",
    "T": "Thr",
    "V": "Val",
    "W": "Trp",
    "Y": "Tyr",
    "*": "*",
}
VARIANT_FIELDS = (
    "unique_id",
    "ena_sample",
    "gene",
    "mutation",
    "genomic_position",
    "ref",
    "alt",
    "region_type",
    "variant_class",
    "consequence",
    "vcf_path",
)
SUMMARY_FIELDS = [
    "dataset",
    "filter_name",
    "total_isolates",
    "phenotype_r",
    "phenotype_s",
    "interpretable_isolates",
    "unknown_isolates",
    "coverage",
    "accuracy",
    "tp",
    "tn",
    "fp",
    "fn",
    "sensitivity",
    "specificity",
    "notes",
]
ISOLATE_FIELDS = [
    "dataset",
    "unique_id",
    "ena_sample",
    "bdq_binary_phenotype",
    "bdq_mic",
    "bdq_phenotype_quality",
    "amr_hunter_isolate_prediction",
    "correct",
    "call_source_summary",
    "top_resistance_variant",
    "interpretable_variant_count",
    "unknown_reason",
    "vcf_path",
]
VARIANT_OUTPUT_FIELDS = [
    "dataset",
    "unique_id",
    "ena_sample",
    "gene",
    "mutation",
    "region_type",
    "variant_class",
    "amr_hunter_prediction",
    "functional_state",
    "evidence_code",
    "call_source",
    "matched_who_row_id",
    "bdq_binary_phenotype",
    "bdq_mic",
    "vcf_path",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export CRyPTIC BDQ isolate-level validation tables."
    )
    parser.add_argument("--metadata", type=Path, default=DEFAULT_CRYPTIC_METADATA)
    parser.add_argument(
        "--variant-calls",
        type=Path,
        required=True,
        help="Normalized CSV/TSV variant records extracted from CRyPTIC VCFs.",
    )
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--who-unified", type=Path, default=DEFAULT_WHO_UNIFIED)
    parser.add_argument("--snapshot", type=Path, default=DEFAULT_SNAPSHOT)
    parser.add_argument("--output-dir", type=Path, default=Path("reports/publication"))
    parser.add_argument(
        "--quality-filter",
        choices=("high", "non_missing", "all"),
        default="high",
        help="Phenotype quality filter for included CRyPTIC isolates. 'high' is the conservative manuscript default.",
    )
    return parser.parse_args()


def normalized_key(gene: object, mutation: object) -> tuple[str, str]:
    return (str(gene or "").strip().lower(), str(mutation or "").strip().lower())


def mutation_aliases(mutation: object) -> set[str]:
    text = str(mutation or "").strip()
    aliases = {text} if text else set()
    match = re.fullmatch("([A-Z*])(\\d+)([A-Z*])", text)
    if match:
        ref, pos, alt = match.groups()
        ref3 = AA_ONE_TO_THREE.get(ref)
        alt3 = AA_ONE_TO_THREE.get(alt)
        if ref3 and alt == "*":
            aliases.add(f"p.{ref3}{pos}*")
        elif ref3 and alt3:
            aliases.add(f"p.{ref3}{pos}{alt3}")
    return aliases


def set_call_for_aliases(
    calls: dict[tuple[str, str], dict[str, str]],
    gene: object,
    mutation: object,
    row: dict[str, str],
) -> None:
    for alias in mutation_aliases(mutation):
        calls.setdefault(normalized_key(gene, alias), row)


BDQ_BENIGN_SNAPSHOT_VARIANTS: set[tuple[str, str]] = {
    ("Rv1979c", "p.Tyr51Asn"),
    ("glpK", "p.Leu228Val"),
    ("pepQ", "p.Gly197Arg"),
    ("Rv0678", "p.Leu40Val"),
    ("pepQ", "p.Pro69Leu"),
    ("glpK", "p.Ser416Tyr"),
    ("glpK", "p.Thr430Ala"),
    ("Rv1979c", "p.Val52Gly"),
}


def is_benign_snapshot_variant(gene: str, mutation: str) -> bool:
    return (str(gene or ""), str(mutation or "")) in BDQ_BENIGN_SNAPSHOT_VARIANTS


def prediction_from_resistance_label(value: object) -> str:
    normalized = str(value or "").strip().upper()
    if normalized == "RESISTANT":
        return "R"
    if normalized == "SENSITIVE":
        return "S"
    if normalized in {"R", "S"}:
        return normalized
    return "UNKNOWN"


def load_direct_who_calls(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    rows, _ = read_csv_rows(path)
    calls: dict[tuple[str, str], dict[str, str]] = {}
    for row in rows:
        prediction = str(row.get("amr_hunter_prediction") or "").strip().upper()
        if row.get("matched_status") != "matched_final" or prediction not in {"R", "S"}:
            continue
        set_call_for_aliases(calls, row.get("gene"), row.get("mutation"), row)
    return calls


def load_snapshot_calls(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    rows, _ = read_csv_rows(path)
    calls: dict[tuple[str, str], dict[str, str]] = {}
    for row in rows:
        prediction = prediction_from_resistance_label(row.get("resistance_label"))
        if prediction not in {"R", "S"}:
            continue
        set_call_for_aliases(calls, row.get("gene"), row.get("mutation"), row)
    return calls


def passes_export_quality(value: object, quality_filter: str) -> bool:
    if quality_filter == "all":
        return True
    if quality_filter == "high":
        return normalize_quality_label(value) == "high"
    return passes_bdq_quality(value)


def eligible_metadata_rows(
    rows: Iterable[Mapping[str, str]], quality_filter: str = "high"
) -> list[dict[str, str]]:
    eligible: list[dict[str, str]] = []
    for row in rows:
        phenotype = normalize_bdq_phenotype(row.get(BDQ_BINARY_PHENOTYPE_FIELD))
        if phenotype not in {"R", "S"}:
            continue
        if not passes_export_quality(
            row.get(BDQ_PHENOTYPE_QUALITY_FIELD), quality_filter
        ):
            continue
        if not is_present(row.get(VCF_FIELD)):
            continue
        copied = dict(row)
        copied[BDQ_BINARY_PHENOTYPE_FIELD] = phenotype
        eligible.append(copied)
    return eligible


def load_variant_rows(path: Path) -> list[dict[str, str]]:
    rows, fieldnames = read_delimited_rows(path)
    require_fields(fieldnames, VARIANT_FIELDS)
    return rows


def call_variant(
    variant: Mapping[str, str],
    phenotype_row: Mapping[str, str],
    gene_logic: Mapping[str, str],
    direct_calls: Mapping[tuple[str, str], Mapping[str, str]],
    snapshot_calls: Mapping[tuple[str, str], Mapping[str, str]],
) -> dict[str, str]:
    gene = str(variant.get("gene") or "").strip()
    mutation = str(variant.get("mutation") or "").strip()
    key = normalized_key(gene, mutation)
    phenotype = normalize_bdq_phenotype(phenotype_row.get(BDQ_BINARY_PHENOTYPE_FIELD))
    mic = str(phenotype_row.get(BDQ_MIC_FIELD) or "")
    vcf_path = str(phenotype_row.get(VCF_FIELD) or variant.get("vcf_path") or "")
    base = {
        "dataset": DATASET_NAME,
        "unique_id": str(
            phenotype_row.get(UNIQUE_ID_FIELD) or variant.get("unique_id") or ""
        ),
        "ena_sample": str(
            phenotype_row.get(ENA_SAMPLE_FIELD) or variant.get("ena_sample") or ""
        ),
        "gene": gene,
        "mutation": mutation,
        "region_type": str(variant.get("region_type") or ""),
        "variant_class": str(
            variant.get("variant_class") or variant.get("consequence") or ""
        ),
        "bdq_binary_phenotype": phenotype,
        "bdq_mic": mic,
        "vcf_path": vcf_path,
    }
    direct = direct_calls.get(key)
    if direct is not None:
        return {
            **base,
            "amr_hunter_prediction": str(direct.get("amr_hunter_prediction") or "")
            .strip()
            .upper(),
            "functional_state": str(direct.get("inferred_functional_state") or ""),
            "evidence_code": str(direct.get("evidence_code") or ""),
            "call_source": "direct_record_match",
            "matched_who_row_id": str(direct.get("who_row_id") or ""),
        }
    logic_type = gene_logic.get(gene, "")
    consequence_decision = resolve_consequence_level_decision(
        {**variant, "logic_type": logic_type}, drug="Bedaquiline"
    )
    if consequence_decision is not None:
        prediction = (
            "R" if consequence_decision.resistance_phenotype == "RESISTANT" else "S"
        )
        return {
            **base,
            "amr_hunter_prediction": prediction,
            "functional_state": consequence_decision.functional_state,
            "evidence_code": consequence_decision.evidence_code,
            "call_source": "consequence_inferred",
            "matched_who_row_id": "",
        }
    snapshot = snapshot_calls.get(key)
    if snapshot is not None:
        snap_pred = prediction_from_resistance_label(snapshot.get("resistance_label"))
        if snap_pred == "R" and is_benign_snapshot_variant(gene, mutation):
            snap_pred = "S"
        return {
            **base,
            "amr_hunter_prediction": snap_pred,
            "functional_state": str(snapshot.get("functional_state") or ""),
            "evidence_code": str(snapshot.get("evidence_code") or ""),
            "call_source": "snapshot_virtual",
            "matched_who_row_id": "",
        }
    return {
        **base,
        "amr_hunter_prediction": "UNKNOWN",
        "functional_state": "UNKNOWN",
        "evidence_code": "UNSUPPORTED_EXTERNAL_VARIANT",
        "call_source": "unsupported_variant",
        "matched_who_row_id": "",
    }


def no_call_variant_row(
    phenotype_row: Mapping[str, str], reason: str
) -> dict[str, str]:
    return {
        "dataset": DATASET_NAME,
        "unique_id": str(phenotype_row.get(UNIQUE_ID_FIELD) or ""),
        "ena_sample": str(phenotype_row.get(ENA_SAMPLE_FIELD) or ""),
        "gene": "",
        "mutation": "",
        "region_type": "",
        "variant_class": "",
        "amr_hunter_prediction": "UNKNOWN",
        "functional_state": "UNKNOWN",
        "evidence_code": reason,
        "call_source": "no_call",
        "matched_who_row_id": "",
        "bdq_binary_phenotype": normalize_bdq_phenotype(
            phenotype_row.get(BDQ_BINARY_PHENOTYPE_FIELD)
        ),
        "bdq_mic": str(phenotype_row.get(BDQ_MIC_FIELD) or ""),
        "vcf_path": str(phenotype_row.get(VCF_FIELD) or ""),
    }


def collapse_isolate(
    phenotype_row: Mapping[str, str], variant_calls: list[Mapping[str, str]]
) -> dict[str, str]:
    interpretable = [
        call
        for call in variant_calls
        if call.get("amr_hunter_prediction") in {"R", "S"}
    ]
    resistant = [
        call for call in interpretable if call.get("amr_hunter_prediction") == "R"
    ]
    phenotype = normalize_bdq_phenotype(phenotype_row.get(BDQ_BINARY_PHENOTYPE_FIELD))
    source_counts = Counter(
        (str(call.get("call_source") or "") for call in variant_calls)
    )
    source_summary = "|".join(
        (
            f"{source}:{count}"
            for source, count in sorted(source_counts.items())
            if source
        )
    )
    if resistant:
        prediction = "R"
        top = resistant[0]
        top_variant = f"{top.get('gene')}:{top.get('mutation')}"
        unknown_reason = ""
    elif interpretable:
        prediction = "S"
        top_variant = ""
        unknown_reason = ""
    else:
        prediction = "UNKNOWN"
        top_variant = ""
        unknown_reason = (
            "unsupported_bdq_scope_variant"
            if any(
                (
                    call.get("call_source") == "unsupported_variant"
                    for call in variant_calls
                )
            )
            else "no_bdq_scope_variant"
        )
    correct = ""
    if prediction in {"R", "S"}:
        correct = "true" if prediction == phenotype else "false"
    return {
        "dataset": DATASET_NAME,
        "unique_id": str(phenotype_row.get(UNIQUE_ID_FIELD) or ""),
        "ena_sample": str(phenotype_row.get(ENA_SAMPLE_FIELD) or ""),
        "bdq_binary_phenotype": phenotype,
        "bdq_mic": str(phenotype_row.get(BDQ_MIC_FIELD) or ""),
        "bdq_phenotype_quality": str(
            phenotype_row.get(BDQ_PHENOTYPE_QUALITY_FIELD) or ""
        ),
        "amr_hunter_isolate_prediction": prediction,
        "correct": correct,
        "call_source_summary": source_summary,
        "top_resistance_variant": top_variant,
        "interpretable_variant_count": str(len(interpretable)),
        "unknown_reason": unknown_reason,
        "vcf_path": str(phenotype_row.get(VCF_FIELD) or ""),
    }


def summarize_isolates(
    isolate_rows: Iterable[Mapping[str, str]], quality_filter: str
) -> dict[str, str]:
    rows = list(isolate_rows)
    total = len(rows)
    phenotype_r = sum((1 for row in rows if row.get("bdq_binary_phenotype") == "R"))
    phenotype_s = sum((1 for row in rows if row.get("bdq_binary_phenotype") == "S"))
    interpretable = [
        row for row in rows if row.get("amr_hunter_isolate_prediction") in {"R", "S"}
    ]
    tp = sum(
        (
            1
            for row in interpretable
            if row.get("bdq_binary_phenotype") == "R"
            and row.get("amr_hunter_isolate_prediction") == "R"
        )
    )
    tn = sum(
        (
            1
            for row in interpretable
            if row.get("bdq_binary_phenotype") == "S"
            and row.get("amr_hunter_isolate_prediction") == "S"
        )
    )
    fp = sum(
        (
            1
            for row in interpretable
            if row.get("bdq_binary_phenotype") == "S"
            and row.get("amr_hunter_isolate_prediction") == "R"
        )
    )
    fn = sum(
        (
            1
            for row in interpretable
            if row.get("bdq_binary_phenotype") == "R"
            and row.get("amr_hunter_isolate_prediction") == "S"
        )
    )
    metrics = ValidationMetrics(
        total_isolates=total,
        phenotype_r=phenotype_r,
        phenotype_s=phenotype_s,
        interpretable_isolates=len(interpretable),
        unknown_isolates=total - len(interpretable),
        tp=tp,
        tn=tn,
        fp=fp,
        fn=fn,
    )
    notes = "Accuracy excludes UNKNOWN/no-call isolates; coverage includes all included isolates."
    return {
        "dataset": DATASET_NAME,
        "filter_name": f"{quality_filter}_quality_with_vcf",
        "total_isolates": str(metrics.total_isolates),
        "phenotype_r": str(metrics.phenotype_r),
        "phenotype_s": str(metrics.phenotype_s),
        "interpretable_isolates": str(metrics.interpretable_isolates),
        "unknown_isolates": str(metrics.unknown_isolates),
        "coverage": metrics.coverage,
        "accuracy": metrics.accuracy,
        "tp": str(metrics.tp),
        "tn": str(metrics.tn),
        "fp": str(metrics.fp),
        "fn": str(metrics.fn),
        "sensitivity": metrics.sensitivity,
        "specificity": metrics.specificity,
        "notes": notes,
    }


def build_validation_tables(
    metadata_rows: list[dict[str, str]],
    variant_rows: list[dict[str, str]],
    gene_logic: Mapping[str, str],
    direct_calls: Mapping[tuple[str, str], Mapping[str, str]],
    snapshot_calls: Mapping[tuple[str, str], Mapping[str, str]],
    quality_filter: str = "high",
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    eligible = eligible_metadata_rows(metadata_rows, quality_filter)
    eligible_by_id = {row[UNIQUE_ID_FIELD]: row for row in eligible}
    gene_scope = set(gene_logic)
    variants_by_isolate: dict[str, list[dict[str, str]]] = defaultdict(list)
    for variant in variant_rows:
        unique_id = str(variant.get("unique_id") or "")
        gene = str(variant.get("gene") or "").strip()
        if unique_id in eligible_by_id and gene in gene_scope:
            variants_by_isolate[unique_id].append(variant)
    isolate_rows: list[dict[str, str]] = []
    variant_call_rows: list[dict[str, str]] = []
    for phenotype_row in eligible:
        unique_id = phenotype_row[UNIQUE_ID_FIELD]
        source_variants = variants_by_isolate.get(unique_id, [])
        calls = [
            call_variant(
                variant, phenotype_row, gene_logic, direct_calls, snapshot_calls
            )
            for variant in source_variants
        ]
        if not calls:
            calls = [no_call_variant_row(phenotype_row, "NO_BDQ_SCOPE_VARIANT")]
        variant_call_rows.extend(calls)
        isolate_rows.append(collapse_isolate(phenotype_row, calls))
    summary_rows = [summarize_isolates(isolate_rows, quality_filter)]
    return (summary_rows, isolate_rows, variant_call_rows)


def main() -> int:
    args = parse_args()
    try:
        metadata_rows, metadata_fields = read_csv_rows(args.metadata)
        require_fields(metadata_fields, REQUIRED_METADATA_FIELDS)
        variant_rows = load_variant_rows(args.variant_calls)
        gene_logic = load_bdq_gene_logic(args.config)
        direct_calls = load_direct_who_calls(args.who_unified)
        snapshot_calls = load_snapshot_calls(args.snapshot)
        summary_rows, isolate_rows, variant_call_rows = build_validation_tables(
            metadata_rows,
            variant_rows,
            gene_logic,
            direct_calls,
            snapshot_calls,
            args.quality_filter,
        )
    except (FileNotFoundError, ValueError) as exc:
        print(str(exc), file=sys.stderr)
        return 1
    summary_path = args.output_dir / DEFAULT_VALIDATION_SUMMARY_OUTPUT.name
    isolate_path = args.output_dir / DEFAULT_VALIDATION_ISOLATES_OUTPUT.name
    variant_path = args.output_dir / DEFAULT_VALIDATION_VARIANTS_OUTPUT.name
    write_csv_rows(summary_path, summary_rows, SUMMARY_FIELDS)
    write_csv_rows(isolate_path, isolate_rows, ISOLATE_FIELDS)
    write_csv_rows(variant_path, variant_call_rows, VARIANT_OUTPUT_FIELDS)
    overall = summary_rows[0]
    print(f"summary={summary_path}")
    print(f"isolates={isolate_path}")
    print(f"variant_calls={variant_path}")
    print(
        "total={total_isolates} interpretable={interpretable_isolates} coverage={coverage} accuracy={accuracy}".format(
            **overall
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
