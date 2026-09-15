from __future__ import annotations
import csv
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
import re
from typing import Iterable, Mapping

DATASET_NAME = "CRyPTIC"
DEFAULT_CRYPTIC_METADATA = Path(
    "data/external/cryptic/CRyPTIC_reuse_table_20240917.csv"
)
DEFAULT_METADATA_INSPECTION_OUTPUT = Path(
    "reports/publication/external_bdq_phenotype_metadata_inspection.csv"
)
DEFAULT_VALIDATION_SUMMARY_OUTPUT = Path(
    "reports/publication/external_bdq_phenotype_validation_summary.csv"
)
DEFAULT_VALIDATION_ISOLATES_OUTPUT = Path(
    "reports/publication/external_bdq_phenotype_validation_isolates.csv"
)
DEFAULT_VALIDATION_VARIANTS_OUTPUT = Path(
    "reports/publication/external_bdq_phenotype_validation_variant_calls.csv"
)
UNIQUE_ID_FIELD = "UNIQUEID"
ENA_SAMPLE_FIELD = "ENA_SAMPLE"
VCF_FIELD = "VCF"
BDQ_BINARY_PHENOTYPE_FIELD = "BDQ_BINARY_PHENOTYPE"
BDQ_MIC_FIELD = "BDQ_MIC"
BDQ_PHENOTYPE_QUALITY_FIELD = "BDQ_PHENOTYPE_QUALITY"
REQUIRED_METADATA_FIELDS = (
    UNIQUE_ID_FIELD,
    BDQ_BINARY_PHENOTYPE_FIELD,
    BDQ_MIC_FIELD,
    BDQ_PHENOTYPE_QUALITY_FIELD,
    VCF_FIELD,
)


@dataclass(frozen=True)
class ValidationMetrics:
    total_isolates: int
    phenotype_r: int
    phenotype_s: int
    interpretable_isolates: int
    unknown_isolates: int
    tp: int
    tn: int
    fp: int
    fn: int

    @property
    def coverage(self) -> str:
        return format_ratio(self.interpretable_isolates, self.total_isolates)

    @property
    def accuracy(self) -> str:
        return format_ratio(self.tp + self.tn, self.interpretable_isolates)

    @property
    def sensitivity(self) -> str:
        return format_ratio(self.tp, self.tp + self.fn)

    @property
    def specificity(self) -> str:
        return format_ratio(self.tn, self.tn + self.fp)


def is_present(value: object) -> bool:
    normalized = str(value or "").strip()
    return bool(normalized) and normalized.upper() not in {
        "NA",
        "N/A",
        "NULL",
        "NONE",
        ".",
    }


def normalize_bdq_phenotype(value: object) -> str:
    normalized = str(value or "").strip().upper()
    if normalized in {"R", "RESISTANT", "RESISTANCE"}:
        return "R"
    if normalized in {"S", "SENSITIVE", "SUSCEPTIBLE", "SUSCEPTIBILITY"}:
        return "S"
    return ""


def normalize_quality_label(value: object) -> str:
    return str(value or "").strip().lower().replace(" ", "_") or "missing"


def passes_bdq_quality(value: object) -> bool:
    label = normalize_quality_label(value)
    if label == "missing":
        return False
    failing_tokens = ("fail", "failed", "poor", "invalid", "low_quality")
    return not any((token in label for token in failing_tokens))


def require_fields(fieldnames: Iterable[str] | None, required: Iterable[str]) -> None:
    present = set(fieldnames or [])
    missing = [field for field in required if field not in present]
    if missing:
        raise ValueError(
            "Missing required CRyPTIC metadata field(s): " + ", ".join(missing)
        )


def read_csv_rows(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        return (rows, list(reader.fieldnames or []))


def read_delimited_rows(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    delimiter = "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        rows = list(reader)
        return (rows, list(reader.fieldnames or []))


def format_ratio(numerator: int, denominator: int) -> str:
    if denominator == 0:
        return ""
    return f"{numerator / denominator:.4f}"


def load_bdq_gene_logic(config_path: Path) -> dict[str, str]:
    try:
        import yaml
    except ModuleNotFoundError:
        return load_bdq_gene_logic_without_yaml(config_path)
    with config_path.open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}
    logic: dict[str, str] = {}
    for gene_cfg in config.get("target_genes", []) or []:
        gene = str(gene_cfg.get("gene") or "").strip()
        logic_type = str(gene_cfg.get("logic_type") or "").strip().upper()
        if gene and logic_type:
            logic[gene] = logic_type
    return logic


def load_bdq_gene_logic_without_yaml(config_path: Path) -> dict[str, str]:
    logic: dict[str, str] = {}
    in_target_genes = False
    current_gene = ""
    current_logic = ""
    key_pattern = re.compile("^([A-Za-z_][\\w-]*):\\s*(.*?)\\s*(?:#.*)?$")

    def flush_current() -> None:
        if current_gene and current_logic:
            logic[current_gene] = current_logic.upper()

    with config_path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            stripped = raw_line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if not raw_line.startswith((" ", "\t", "-")) and stripped.endswith(":"):
                if in_target_genes:
                    flush_current()
                    current_gene = ""
                    current_logic = ""
                in_target_genes = stripped == "target_genes:"
                continue
            if not in_target_genes:
                continue
            entry = stripped[2:].strip() if stripped.startswith("- ") else stripped
            if stripped.startswith("- "):
                flush_current()
                current_gene = ""
                current_logic = ""
            match = key_pattern.match(entry)
            if match is None:
                continue
            key, value = match.groups()
            value = value.strip().strip("'\"")
            if key == "gene":
                current_gene = value
            elif key == "logic_type":
                current_logic = value.upper()
    if in_target_genes:
        flush_current()
    return logic


def summarize_metadata_rows(rows: Iterable[Mapping[str, str]]) -> list[dict[str, str]]:
    total = 0
    phenotype_counts: Counter[str] = Counter()
    quality_counts: Counter[str] = Counter()
    mic_available = 0
    vcf_available = 0
    quality_passing = Counter()
    quality_passing_with_vcf = Counter()
    high_quality = Counter()
    high_quality_with_vcf = Counter()
    for row in rows:
        total += 1
        phenotype = normalize_bdq_phenotype(row.get(BDQ_BINARY_PHENOTYPE_FIELD))
        phenotype_counts[phenotype or "missing_or_ambiguous"] += 1
        quality = normalize_quality_label(row.get(BDQ_PHENOTYPE_QUALITY_FIELD))
        quality_counts[quality] += 1
        has_mic = is_present(row.get(BDQ_MIC_FIELD))
        has_vcf = is_present(row.get(VCF_FIELD))
        if has_mic:
            mic_available += 1
        if has_vcf:
            vcf_available += 1
        if phenotype in {"R", "S"} and passes_bdq_quality(
            row.get(BDQ_PHENOTYPE_QUALITY_FIELD)
        ):
            quality_passing[phenotype] += 1
            if has_vcf:
                quality_passing_with_vcf[phenotype] += 1
        if phenotype in {"R", "S"} and quality == "high":
            high_quality[phenotype] += 1
            if has_vcf:
                high_quality_with_vcf[phenotype] += 1
    rows_out = [
        metric_row("total_rows", total),
        metric_row("bdq_binary_phenotype_r", phenotype_counts["R"]),
        metric_row("bdq_binary_phenotype_s", phenotype_counts["S"]),
        metric_row(
            "bdq_binary_phenotype_missing", phenotype_counts["missing_or_ambiguous"]
        ),
        metric_row("bdq_mic_available", mic_available),
        metric_row("bdq_mic_missing", total - mic_available),
        metric_row("vcf_available", vcf_available),
        metric_row("vcf_missing", total - vcf_available),
        metric_row("quality_passing_r", quality_passing["R"]),
        metric_row("quality_passing_s", quality_passing["S"]),
        metric_row("quality_passing_with_vcf_r", quality_passing_with_vcf["R"]),
        metric_row("quality_passing_with_vcf_s", quality_passing_with_vcf["S"]),
        metric_row("high_quality_r", high_quality["R"]),
        metric_row("high_quality_s", high_quality["S"]),
        metric_row("high_quality_with_vcf_r", high_quality_with_vcf["R"]),
        metric_row("high_quality_with_vcf_s", high_quality_with_vcf["S"]),
    ]
    for quality, count in sorted(quality_counts.items()):
        rows_out.append(metric_row(f"bdq_quality_{quality}", count))
    feasible = quality_passing_with_vcf["R"] > 0 and quality_passing_with_vcf["S"] > 0
    note = (
        "At least one quality-passing R and S isolate with VCF is available."
        if feasible
        else "Insufficient quality-passing R/S isolates with VCF for balanced validation."
    )
    rows_out.append(
        metric_row(
            "validation_feasibility", "feasible" if feasible else "insufficient", note
        )
    )
    return rows_out


def metric_row(metric: str, value: object, notes: str = "") -> dict[str, str]:
    return {
        "dataset": DATASET_NAME,
        "metric": metric,
        "value": str(value),
        "notes": notes,
    }


def write_metadata_summary(path: Path, rows: Iterable[Mapping[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["dataset", "metric", "value", "notes"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def write_csv_rows(
    path: Path, rows: Iterable[Mapping[str, object]], fieldnames: list[str]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})
