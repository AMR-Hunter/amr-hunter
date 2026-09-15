from __future__ import annotations
import argparse
import csv
from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor, as_completed
import gzip
from pathlib import Path
import re
import sys
import time
from typing import Iterable, Mapping
from urllib.request import urlretrieve

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
from scripts.cryptic_bdq_validation_common import (
    BDQ_BINARY_PHENOTYPE_FIELD,
    BDQ_PHENOTYPE_QUALITY_FIELD,
    DEFAULT_CRYPTIC_METADATA,
    ENA_SAMPLE_FIELD,
    REQUIRED_METADATA_FIELDS,
    UNIQUE_ID_FIELD,
    VCF_FIELD,
    is_present,
    load_bdq_gene_logic,
    normalize_bdq_phenotype,
    normalize_quality_label,
    passes_bdq_quality,
    read_csv_rows,
    require_fields,
)

DEFAULT_CONFIG = Path("config/config.yaml")
DEFAULT_REFERENCE_GBK = Path("data/reference/mtb_h37rv.gbk")
DEFAULT_OUTPUT = Path("data/external/cryptic/bdq_scope_variant_calls.tsv")
DEFAULT_VCF_ROOT = Path("data/external/cryptic")
DEFAULT_VCF_BASE_URL = "https://ftp.ebi.ac.uk/pub/databases/cryptic/release_june2022/"
OUTPUT_FIELDS = [
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
]
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
    "X": "Xaa",
}
CODON_TABLE_11 = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "W",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


@dataclass(frozen=True)
class GeneFeature:
    gene: str
    start: int
    end: int
    strand: int
    sequence: str
    upstream_bp: int = 120


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract normalized BDQ-scope variants from CRyPTIC VCFs."
    )
    parser.add_argument("--metadata", type=Path, default=DEFAULT_CRYPTIC_METADATA)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--reference-gbk", type=Path, default=DEFAULT_REFERENCE_GBK)
    parser.add_argument("--vcf-root", type=Path, default=DEFAULT_VCF_ROOT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--vcf-base-url", default=DEFAULT_VCF_BASE_URL)
    parser.add_argument(
        "--quality-filter", choices=("high", "non_missing", "all"), default="high"
    )
    parser.add_argument("--download-missing", action="store_true")
    parser.add_argument("--max-isolates", type=int, default=0)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument(
        "--genes",
        default="",
        help="Optional comma-separated gene subset for tests or targeted extraction.",
    )
    return parser.parse_args()


def sanitize_sequence(value: object) -> str:
    return str(value or "").upper().replace("U", "T")


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(table)[::-1].upper()


def translate_codon(codon: str) -> str:
    return CODON_TABLE_11.get(sanitize_sequence(codon), "X")


def aa3(value: str) -> str:
    return AA_ONE_TO_THREE.get(value, "Xaa")


def normalize_genbank_location(value: str) -> tuple[int, int, int]:
    text = value.replace("<", "").replace(">", "")
    strand = -1 if "complement" in text else 1
    numbers = [int(item) for item in re.findall("\\d+", text)]
    if len(numbers) < 2:
        raise ValueError(f"Unsupported GenBank location: {value}")
    return (min(numbers), max(numbers), strand)


def parse_genbank(path: Path, target_genes: Iterable[str]) -> dict[str, GeneFeature]:
    target = set(target_genes)
    sequence_chunks: list[str] = []
    feature_blocks: list[list[str]] = []
    current_block: list[str] | None = None
    in_features = False
    in_origin = False
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("FEATURES"):
                in_features = True
                continue
            if line.startswith("ORIGIN"):
                if current_block is not None:
                    feature_blocks.append(current_block)
                    current_block = None
                in_features = False
                in_origin = True
                continue
            if in_origin:
                if line.startswith("//"):
                    break
                sequence_chunks.append("".join(re.findall("[A-Za-z]+", line)))
                continue
            if not in_features:
                continue
            feature_type = line[5:21].strip() if len(line) >= 21 else ""
            if feature_type:
                if current_block is not None:
                    feature_blocks.append(current_block)
                current_block = [line.rstrip("\n")]
            elif current_block is not None:
                current_block.append(line.rstrip("\n"))
    genome = sanitize_sequence("".join(sequence_chunks))
    if not genome:
        raise ValueError(f"No ORIGIN sequence found in {path}")
    features: dict[str, GeneFeature] = {}
    for block in feature_blocks:
        first = block[0]
        feature_type = first[5:21].strip() if len(first) >= 21 else ""
        if feature_type != "CDS":
            continue
        qualifiers: dict[str, list[str]] = {}
        for line in block[1:]:
            stripped = line.strip()
            if not stripped.startswith("/"):
                continue
            key, _, raw_value = stripped[1:].partition("=")
            qualifiers.setdefault(key, []).append(raw_value.strip().strip('"'))
        names = {
            value
            for key in ("gene", "locus_tag", "old_locus_tag")
            for value in qualifiers.get(key, [])
        }
        matched = sorted(names & target)
        if not matched:
            continue
        gene = matched[0]
        start, end, strand = normalize_genbank_location(first[21:].strip())
        genomic_seq = genome[start - 1 : end]
        cds_seq = genomic_seq if strand == 1 else reverse_complement(genomic_seq)
        features[gene] = GeneFeature(
            gene=gene, start=start, end=end, strand=strand, sequence=cds_seq
        )
    missing = target - set(features)
    if missing:
        raise ValueError(
            f"Missing GenBank CDS features for configured genes: {sorted(missing)}"
        )
    return features


def passes_export_quality(value: object, quality_filter: str) -> bool:
    if quality_filter == "all":
        return True
    if quality_filter == "high":
        return normalize_quality_label(value) == "high"
    return passes_bdq_quality(value)


def eligible_metadata_rows(
    rows: Iterable[Mapping[str, str]], quality_filter: str
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


def transcript_ref_alt(feature: GeneFeature, ref: str, alt: str) -> tuple[str, str]:
    if feature.strand == 1:
        return (sanitize_sequence(ref), sanitize_sequence(alt))
    return (reverse_complement(ref), reverse_complement(alt))


def cds_position(feature: GeneFeature, genomic_pos: int, ref: str) -> int:
    if feature.strand == 1:
        return genomic_pos - feature.start + 1
    return feature.end - (genomic_pos + len(ref) - 1) + 1


def annotate_cds_variant(
    feature: GeneFeature, genomic_pos: int, ref: str, alt: str
) -> tuple[str, str, str] | None:
    ref_t, alt_t = transcript_ref_alt(feature, ref, alt)
    pos = cds_position(feature, genomic_pos, ref)
    if pos < 1 or pos > len(feature.sequence):
        return None
    observed = feature.sequence[pos - 1 : pos - 1 + len(ref_t)]
    if observed != ref_t:
        return None
    if len(ref_t) == 1 and len(alt_t) == 1:
        codon_start = (pos - 1) // 3 * 3
        codon_offset = (pos - 1) % 3
        wt_codon = feature.sequence[codon_start : codon_start + 3]
        mt_codon = wt_codon[:codon_offset] + alt_t + wt_codon[codon_offset + 1 :]
        wt_aa = translate_codon(wt_codon)
        mt_aa = translate_codon(mt_codon)
        aa_pos = codon_start // 3 + 1
        if aa_pos == 1 and wt_aa == "M" and (mt_aa != "M"):
            return ("p.Met1?", "start_lost", "start_lost")
        if wt_aa == mt_aa:
            return (f"c.{pos}{ref_t}>{alt_t}", "snv", "synonymous_variant")
        if mt_aa == "*":
            return (f"p.{aa3(wt_aa)}{aa_pos}*", "nonsense", "stop_gained")
        return (
            f"p.{aa3(wt_aa)}{aa_pos}{aa3(mt_aa)}",
            "protein_substitution",
            "missense_variant",
        )
    aa_pos = (pos - 1) // 3 + 1
    wt_codon_start = (pos - 1) // 3 * 3
    wt_aa = translate_codon(feature.sequence[wt_codon_start : wt_codon_start + 3])
    if (len(alt_t) - len(ref_t)) % 3 != 0:
        return (f"p.{aa3(wt_aa)}{aa_pos}fs", "frameshift", "frameshift_variant")
    return (f"c.{pos}{ref_t}>{alt_t}", "inframe_indel", "inframe_indel")


def regulatory_label(
    feature: GeneFeature, genomic_pos: int, ref: str, alt: str
) -> str | None:
    ref_t, alt_t = transcript_ref_alt(feature, ref, alt)
    if feature.strand == 1 and genomic_pos < feature.start:
        offset = feature.start - genomic_pos
    elif feature.strand == -1 and genomic_pos > feature.end:
        offset = genomic_pos - feature.end
    else:
        return None
    return f"c.-{offset}{ref_t}>{alt_t}"


def in_regulatory_window(feature: GeneFeature, genomic_pos: int) -> bool:
    if feature.start <= genomic_pos <= feature.end:
        return False
    if feature.strand == 1:
        return feature.start - feature.upstream_bp <= genomic_pos < feature.start
    return feature.end < genomic_pos <= feature.end + feature.upstream_bp


def annotate_variant(
    features: Mapping[str, GeneFeature], genomic_pos: int, ref: str, alt: str
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for feature in features.values():
        if feature.start <= genomic_pos <= feature.end:
            annotated = annotate_cds_variant(feature, genomic_pos, ref, alt)
            if annotated is None:
                continue
            mutation, variant_class, consequence = annotated
            rows.append(
                {
                    "gene": feature.gene,
                    "mutation": mutation,
                    "region_type": "CDS",
                    "variant_class": variant_class,
                    "consequence": consequence,
                }
            )
        elif (
            in_regulatory_window(feature, genomic_pos)
            and len(ref) == 1
            and (len(alt) == 1)
        ):
            mutation = regulatory_label(feature, genomic_pos, ref, alt)
            if mutation is None:
                continue
            rows.append(
                {
                    "gene": feature.gene,
                    "mutation": mutation,
                    "region_type": "REGULATORY",
                    "variant_class": "regulatory",
                    "consequence": "upstream_gene_variant",
                }
            )
    return rows


def open_vcf(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open(encoding="utf-8", errors="replace")


def selected_alt_indices(format_value: str, sample_value: str) -> set[int]:
    keys = format_value.split(":")
    values = sample_value.split(":")
    sample = dict(zip(keys, values))
    gt = sample.get("GT", "")
    if not gt or gt in {"0/0", "0|0", "./.", ".|."}:
        return set()
    return {
        int(item) for item in re.split("[/|]", gt) if item.isdigit() and int(item) > 0
    }


def vcf_local_path(vcf_root: Path, relative_path: str) -> Path:
    normalized = relative_path.replace("\\", "/")
    while normalized.startswith("../"):
        normalized = normalized[3:]
    return vcf_root / normalized


def vcf_url(base_url: str, relative_path: str) -> str:
    normalized = relative_path.replace("\\", "/")
    while normalized.startswith("../"):
        normalized = normalized[3:]
    return base_url.rstrip("/") + "/" + normalized


def download_with_retries(url: str, path: Path, attempts: int = 3) -> None:
    last_error: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            urlretrieve(url, path)
            return
        except Exception as exc:
            last_error = exc
            if path.exists():
                path.unlink()
            if attempt < attempts:
                time.sleep(2 * attempt)
    if last_error is not None:
        raise last_error


def ensure_vcf(
    relative_path: str, vcf_root: Path, base_url: str, download_missing: bool
) -> Path:
    path = vcf_local_path(vcf_root, relative_path)
    if path.exists():
        return path
    if not download_missing:
        raise FileNotFoundError(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    download_with_retries(vcf_url(base_url, relative_path), path)
    return path


def extract_vcf_rows(
    phenotype_row: Mapping[str, str],
    features: Mapping[str, GeneFeature],
    vcf_path: Path,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with open_vcf(vcf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            _chrom, pos_text, _id, ref, alts, _qual, filt, _info, fmt, sample = parts[
                :10
            ]
            if filt != "PASS":
                continue
            try:
                genomic_pos = int(pos_text)
            except ValueError:
                continue
            alt_values = alts.split(",")
            for alt_index in selected_alt_indices(fmt, sample):
                if alt_index < 1 or alt_index > len(alt_values):
                    continue
                alt = alt_values[alt_index - 1]
                if alt.startswith("<") or ref.startswith("<"):
                    continue
                for annotation in annotate_variant(features, genomic_pos, ref, alt):
                    rows.append(
                        {
                            "unique_id": phenotype_row.get(UNIQUE_ID_FIELD, ""),
                            "ena_sample": phenotype_row.get(ENA_SAMPLE_FIELD, ""),
                            "genomic_position": str(genomic_pos),
                            "ref": sanitize_sequence(ref),
                            "alt": sanitize_sequence(alt),
                            "vcf_path": str(phenotype_row.get(VCF_FIELD) or vcf_path),
                            **annotation,
                        }
                    )
    return rows


def write_tsv(path: Path, rows: Iterable[Mapping[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in OUTPUT_FIELDS})


def process_isolate(
    row: Mapping[str, str],
    features: Mapping[str, GeneFeature],
    vcf_root: Path,
    vcf_base_url: str,
    download_missing: bool,
) -> tuple[list[dict[str, str]], int, int, str]:
    try:
        path = ensure_vcf(
            str(row.get(VCF_FIELD) or ""), vcf_root, vcf_base_url, download_missing
        )
        return (extract_vcf_rows(row, features, path), 0, 0, "")
    except FileNotFoundError:
        return ([], 1, 0, "")
    except Exception as exc:
        return (
            [],
            0,
            1,
            f"warning: failed to parse VCF for {row.get(UNIQUE_ID_FIELD)}: {exc}",
        )


def main() -> int:
    args = parse_args()
    try:
        metadata_rows, metadata_fields = read_csv_rows(args.metadata)
        require_fields(metadata_fields, REQUIRED_METADATA_FIELDS)
        gene_logic = load_bdq_gene_logic(args.config)
        selected_genes = [
            gene.strip() for gene in str(args.genes or "").split(",") if gene.strip()
        ]
        target_genes = selected_genes or list(gene_logic.keys())
        features = parse_genbank(args.reference_gbk, target_genes)
        eligible = eligible_metadata_rows(metadata_rows, args.quality_filter)
        if args.max_isolates > 0:
            eligible = eligible[: args.max_isolates]
        extracted: list[dict[str, str]] = []
        missing_vcfs = 0
        failed_vcfs = 0
        workers = max(1, int(args.workers or 1))
        if workers == 1:
            for index, row in enumerate(eligible, start=1):
                rows, missing, failed, warning = process_isolate(
                    row,
                    features,
                    args.vcf_root,
                    args.vcf_base_url,
                    args.download_missing,
                )
                extracted.extend(rows)
                missing_vcfs += missing
                failed_vcfs += failed
                if warning:
                    print(warning, file=sys.stderr)
                if index % 500 == 0:
                    print(
                        f"processed={index} extracted_rows={len(extracted)}",
                        file=sys.stderr,
                    )
        else:
            with ThreadPoolExecutor(max_workers=workers) as executor:
                futures = [
                    executor.submit(
                        process_isolate,
                        row,
                        features,
                        args.vcf_root,
                        args.vcf_base_url,
                        args.download_missing,
                    )
                    for row in eligible
                ]
                for index, future in enumerate(as_completed(futures), start=1):
                    rows, missing, failed, warning = future.result()
                    extracted.extend(rows)
                    missing_vcfs += missing
                    failed_vcfs += failed
                    if warning:
                        print(warning, file=sys.stderr)
                    if index % 500 == 0:
                        print(
                            f"processed={index} extracted_rows={len(extracted)}",
                            file=sys.stderr,
                        )
    except (FileNotFoundError, ValueError) as exc:
        print(str(exc), file=sys.stderr)
        return 1
    write_tsv(args.output, extracted)
    print(f"output={args.output}")
    print(f"eligible_isolates={len(eligible)}")
    print(f"variant_rows={len(extracted)}")
    print(f"missing_vcfs={missing_vcfs}")
    print(f"failed_vcfs={failed_vcfs}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
