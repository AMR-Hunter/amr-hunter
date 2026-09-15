from __future__ import annotations
import argparse
import json
import re
import shutil
import sqlite3
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple
import pandas as pd
import yaml
from Bio import SeqIO
from Bio.Seq import Seq

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
from scripts.evaluate_who_catalog import (
    GENERIC_NON_EXECUTABLE_MUTATIONS,
    classify_current_support,
    load_who_catalog,
    normalize_window_specs,
    regulatory_candidate_keys,
    resolve_modeled_genes,
)
from scripts.init_project import extract_target_gene_contexts

DEFAULT_DB = Path("/root/gpufree-data/db/amr_hunter.db")
DEFAULT_XLSX = Path(
    "/root/gpufree-share/amr_hunter_workspace/amr_hunter_data/reference/who/WHO-UCN-TB-2023.5-eng.xlsx"
)
DEFAULT_GBK = Path(
    "/root/gpufree-share/amr_hunter_workspace/amr_hunter_data/reference/mtb_h37rv.gbk"
)
DEFAULT_CONFIG = Path(
    "/root/gpufree-share/amr_hunter_workspace/code/config/config.yaml"
)
DEFAULT_BENCHMARKS_ROOT = Path("/root/gpufree-data/amr_hunter_runtime/benchmarks")
SIMPLE_CDS_MUTATION_RE = re.compile("c\\.(\\d+)([ACGT]+)>([ACGT]+)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create an isolated WHO benchmark DB/config for current-supported truth-labeled rows."
    )
    parser.add_argument(
        "--db", type=Path, default=DEFAULT_DB, help="Source SQLite DB path."
    )
    parser.add_argument(
        "--xlsx", type=Path, default=DEFAULT_XLSX, help="WHO workbook path."
    )
    parser.add_argument(
        "--gbk", type=Path, default=DEFAULT_GBK, help="H37Rv GenBank path."
    )
    parser.add_argument(
        "--config", type=Path, default=DEFAULT_CONFIG, help="Template config.yaml path."
    )
    parser.add_argument(
        "--benchmarks-root",
        type=Path,
        default=DEFAULT_BENCHMARKS_ROOT,
        help="Parent directory under which the benchmark root will be created.",
    )
    parser.add_argument(
        "--benchmark-root",
        type=Path,
        default=None,
        help="Explicit benchmark root. Defaults to benchmarks-root/who_eval_<timestamp>.",
    )
    parser.add_argument(
        "--drug", default="Bedaquiline", help="Drug name in the WHO workbook."
    )
    parser.add_argument(
        "--truth",
        choices=["all", "resistant", "sensitive"],
        default="all",
        help="Restrict truth rows before current-supported filtering.",
    )
    parser.add_argument(
        "--extra-modeled-gene",
        action="append",
        default=[],
        help="Append an optional gene to the default modeled-gene set. Repeatable.",
    )
    return parser.parse_args()


def sanitize_sequence(value: object) -> str:
    return str(value or "").upper().replace("U", "T")


def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def default_benchmark_root(base_root: Path) -> Path:
    timestamp = time.strftime("who_eval_%Y%m%dT%H%M%S")
    return base_root / timestamp


def load_coordinate_rows(xlsx_path: Path) -> Dict[str, List[Dict[str, object]]]:
    coords = pd.read_excel(
        xlsx_path,
        sheet_name="Genomic_coordinates",
        usecols=[
            "variant",
            "position",
            "reference_nucleotide",
            "alternative_nucleotide",
        ],
    )
    coords = coords.dropna(
        subset=["variant", "position", "reference_nucleotide", "alternative_nucleotide"]
    ).copy()
    grouped: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    for row in coords.itertuples(index=False):
        grouped[str(row.variant).strip()].append(
            {
                "position": int(row.position),
                "ref": sanitize_sequence(row.reference_nucleotide),
                "alt": sanitize_sequence(row.alternative_nucleotide),
            }
        )
    return grouped


def load_gene_features(
    gbk_path: Path, genes: Iterable[str]
) -> Dict[str, Dict[str, int]]:
    genes = set(genes)
    features: Dict[str, Dict[str, int]] = {}
    for record in SeqIO.parse(str(gbk_path), "genbank"):
        for feature in record.features:
            if feature.type != "CDS":
                continue
            names = {
                str(value)
                for key in ("gene", "locus_tag", "old_locus_tag")
                for value in feature.qualifiers.get(key, [])
            }
            matched = names & genes
            if not matched:
                continue
            gene_name = sorted(matched)[0]
            features[gene_name] = {
                "start": int(feature.location.start) + 1,
                "end": int(feature.location.end),
                "strand": int(feature.location.strand or 1),
            }
    missing = genes - set(features)
    if missing:
        raise ValueError(f"Missing GenBank features for genes: {sorted(missing)}")
    return features


def normalize_simple_cds_variant(mutation: str) -> Optional[Tuple[int, str, str]]:
    match = SIMPLE_CDS_MUTATION_RE.fullmatch(mutation)
    if not match:
        return None
    return (
        int(match.group(1)),
        sanitize_sequence(match.group(2)),
        sanitize_sequence(match.group(3)),
    )


def convert_genomic_variant_to_cds(
    gene_name: str,
    gene_sequence: str,
    feature: Dict[str, int],
    genomic_pos: int,
    genomic_ref: str,
    genomic_alt: str,
) -> Tuple[int, str, str]:
    start = int(feature["start"])
    end = int(feature["end"])
    strand = int(feature["strand"])
    genomic_ref = sanitize_sequence(genomic_ref)
    genomic_alt = sanitize_sequence(genomic_alt)
    if strand == 1:
        cds_pos = genomic_pos - start + 1
        cds_ref = genomic_ref
        cds_alt = genomic_alt
    else:
        cds_pos = end - (genomic_pos + len(genomic_ref) - 1) + 1
        cds_ref = reverse_complement(genomic_ref)
        cds_alt = reverse_complement(genomic_alt)
    if cds_pos < 1 or cds_pos + len(cds_ref) - 1 > len(gene_sequence):
        raise ValueError(
            f"{gene_name}: CDS edit {cds_pos} {cds_ref}>{cds_alt} falls outside gene sequence length {len(gene_sequence)}"
        )
    observed = gene_sequence[cds_pos - 1 : cds_pos - 1 + len(cds_ref)]
    if cds_ref and observed != cds_ref:
        raise ValueError(
            f"{gene_name}: reference mismatch at CDS {cds_pos}: expected {cds_ref}, observed {observed}"
        )
    return (cds_pos, cds_ref, cds_alt)


def ensure_bedaquiline_id(conn: sqlite3.Connection, drug: str) -> int:
    row = conn.execute("SELECT id FROM drugs WHERE name = ?", (drug,)).fetchone()
    if row is not None:
        return int(row[0])
    conn.execute(
        "INSERT INTO drugs (name, description, created_at, updated_at) VALUES (?, ?, CURRENT_TIMESTAMP, CURRENT_TIMESTAMP)",
        (drug, f"Imported automatically for WHO evaluation benchmark ({drug})"),
    )
    row = conn.execute("SELECT id FROM drugs WHERE name = ?", (drug,)).fetchone()
    if row is None:
        raise RuntimeError(f"Failed to create drug row for {drug}")
    return int(row[0])


def load_benchmark_gene_records(
    conn: sqlite3.Connection,
) -> Dict[str, Dict[str, object]]:
    rows = conn.execute("SELECT id, name, sequence FROM genes").fetchall()
    records: Dict[str, Dict[str, object]] = {}
    for row in rows:
        records[str(row[1])] = {
            "id": int(row[0]),
            "sequence": sanitize_sequence(row[2]),
        }
    return records


def create_benchmark_layout(root: Path) -> Dict[str, Path]:
    paths = {
        "root": root,
        "db_dir": root / "db",
        "results_dir": root / "results",
        "logs_dir": root / "logs",
        "reports_dir": root / "reports",
        "config_path": root / "config.yaml",
        "db_path": root / "db" / "amr_hunter.db",
        "export_csv_path": root / "results" / "amr_results.csv",
        "analysis_snapshot_csv_path": root / "reports" / "analysis_snapshot.csv",
        "log_path": root / "logs" / "pipeline.log",
        "metrics_path": root / "logs" / "batch_metrics.jsonl",
        "manifest_path": root / "reports" / "who_benchmark_manifest.json",
    }
    for key in ("db_dir", "results_dir", "logs_dir", "reports_dir"):
        paths[key].mkdir(parents=True, exist_ok=True)
    return paths


def write_benchmark_config(
    template_path: Path, output_path: Path, paths: Dict[str, Path]
) -> None:
    config = yaml.safe_load(template_path.read_text(encoding="utf-8"))
    paths_cfg = config.setdefault("paths", {})
    paths_cfg["local_db_path"] = str(paths["db_path"])
    paths_cfg["share_db_path"] = str(paths["db_path"])
    pipeline_cfg = config.setdefault("pipeline", {})
    pipeline_cfg["export_csv_path"] = str(paths["export_csv_path"])
    pipeline_cfg["analysis_snapshot_csv_path"] = str(
        paths["analysis_snapshot_csv_path"]
    )
    pipeline_cfg["log_path"] = str(paths["log_path"])
    pipeline_cfg["metrics_path"] = str(paths["metrics_path"])
    output_path.write_text(
        yaml.safe_dump(config, sort_keys=False, allow_unicode=True), encoding="utf-8"
    )


def build_gene_meta_by_name(
    config: Dict[str, Any], target_genes: Iterable[str]
) -> Dict[str, Dict[str, Any]]:
    genes = set(target_genes)
    meta: Dict[str, Dict[str, Any]] = {}
    for item in config.get("target_genes") or []:
        gene = str(item.get("gene") or "").strip()
        if not gene or gene not in genes:
            continue
        meta[gene] = {
            "apply_default_windows": bool(item.get("apply_default_windows", False)),
            "include_upstream": bool(item.get("include_upstream", False)),
            "upstream_bp": int(item.get("upstream_bp", 120) or 120),
            "noncoding_windows": list(item.get("noncoding_windows") or []),
        }
    return meta


def build_regulatory_sequences_by_gene(
    config_path: Path, gbk_path: Path, genes: Iterable[str]
) -> Dict[str, Dict[str, str]]:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    default_windows = list((config.get("noncoding") or {}).get("default_windows") or [])
    intergenic_regions = list(
        (config.get("noncoding") or {}).get("intergenic_regions") or []
    )
    gene_meta_by_name = build_gene_meta_by_name(config, genes)
    contexts = extract_target_gene_contexts(
        gbk_path, set(genes), gene_meta_by_name, default_windows, intergenic_regions
    )
    return {
        gene: {
            str(region.get("name") or "REGULATORY"): str(region.get("sequence") or "")
            for region in context.get("regulatory_regions", [])
            if str(region.get("name") or "").strip()
        }
        for gene, context in contexts.items()
    }


def reset_mutations_table(conn: sqlite3.Connection) -> None:
    conn.execute("DELETE FROM mutations")
    conn.execute("DELETE FROM sqlite_sequence WHERE name = 'mutations'")


def import_variants(
    conn: sqlite3.Connection, rows: Sequence[Tuple[object, ...]]
) -> int:
    conn.executemany(
        "\n        INSERT INTO mutations (\n            gene_id, pos, ref, alt, aa_change, region_type, region_name, source_sequence,\n            evo_model, target_drug_id, status\n        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)\n        ",
        rows,
    )
    result = conn.execute("SELECT COUNT(*) FROM mutations").fetchone()
    return int(result[0]) if result is not None else 0


def select_supported_who_rows(
    xlsx_path: Path,
    drug: str,
    truth_mode: str,
    config_path: Path,
    modeled_genes: Iterable[str],
) -> Tuple[pd.DataFrame, Dict[str, List[Dict[str, Any]]]]:
    gene_window_specs = normalize_window_specs(config_path)
    who = load_who_catalog(xlsx_path, drug, truth_mode, gene_window_specs)
    modeled_gene_set = set(modeled_genes)
    comparable = who[
        who["gene"].isin(modeled_gene_set)
        & who["truth"].notna()
        & who["current_supported"]
    ].copy()
    return (comparable, gene_window_specs)


def build_exact_variants(
    who_rows: pd.DataFrame,
    coordinate_rows: Dict[str, List[Dict[str, object]]],
    gene_records: Dict[str, Dict[str, object]],
    gene_features: Dict[str, Dict[str, int]],
    regulatory_sequences: Dict[str, Dict[str, str]],
    gene_window_specs: Dict[str, List[Dict[str, Any]]],
    drug_id: int,
) -> Tuple[List[Tuple[object, ...]], Dict[str, object]]:
    inserts: List[Tuple[object, ...]] = []
    imported_by_gene: Counter[str] = Counter()
    imported_by_region_type: Counter[str] = Counter()
    who_row_realizations: Counter[str] = Counter()
    generic_rows: List[Dict[str, str]] = []
    unresolved_rows: List[Dict[str, str]] = []
    seen_rows: set[Tuple[object, ...]] = set()
    for _, row in who_rows.iterrows():
        gene = str(row["gene"])
        mutation = str(row["mutation"])
        effect = str(row["effect"] or "")
        gene_record = gene_records.get(gene)
        if gene_record is None:
            unresolved_rows.append(
                {"gene": gene, "mutation": mutation, "reason": "gene_not_modeled_in_db"}
            )
            continue
        gene_sequence = str(gene_record["sequence"])
        if not gene_sequence:
            unresolved_rows.append(
                {"gene": gene, "mutation": mutation, "reason": "missing_gene_sequence"}
            )
            continue
        if mutation in GENERIC_NON_EXECUTABLE_MUTATIONS:
            generic_rows.append(
                {"gene": gene, "mutation": mutation, "reason": "generic_category"}
            )
            continue
        realized = 0
        mapping_errors: List[str] = []
        if effect == "upstream_gene_variant":
            candidate_keys = sorted(
                regulatory_candidate_keys(gene, mutation, gene_window_specs)
            )
            if not candidate_keys:
                unresolved_rows.append(
                    {
                        "gene": gene,
                        "mutation": mutation,
                        "reason": "no_regulatory_candidate_key",
                    }
                )
                continue
            for _, region_name, pos, ref, alt in candidate_keys:
                source_sequence = str(
                    (regulatory_sequences.get(gene) or {}).get(region_name) or ""
                )
                if not source_sequence:
                    mapping_errors.append(f"missing_regulatory_sequence:{region_name}")
                    continue
                if pos < 1 or pos > len(source_sequence):
                    mapping_errors.append(
                        f"regulatory_pos_out_of_bounds:{region_name}:{pos}"
                    )
                    continue
                observed = source_sequence[pos - 1 : pos - 1 + len(ref)]
                if observed != ref:
                    mapping_errors.append(
                        f"regulatory_ref_mismatch:{region_name}:{pos}:{ref}!={observed}"
                    )
                    continue
                insert_row = (
                    int(gene_record["id"]),
                    pos,
                    ref,
                    alt,
                    mutation,
                    "REGULATORY",
                    region_name,
                    source_sequence,
                    "evo2_7b",
                    drug_id,
                    "PENDING",
                )
                if insert_row in seen_rows:
                    continue
                seen_rows.add(insert_row)
                inserts.append(insert_row)
                realized += 1
                imported_by_gene[gene] += 1
                imported_by_region_type["REGULATORY"] += 1
            if realized == 0:
                unresolved_rows.append(
                    {
                        "gene": gene,
                        "mutation": mutation,
                        "reason": (
                            mapping_errors[0]
                            if mapping_errors
                            else "no_regulatory_realization"
                        ),
                    }
                )
                continue
            who_row_realizations[f"{gene}:{mutation}"] = realized
            continue
        exact_variants: List[Tuple[int, str, str]] = []
        simple_cds = normalize_simple_cds_variant(mutation)
        if simple_cds is not None:
            cds_pos, cds_ref, cds_alt = simple_cds
            observed = gene_sequence[cds_pos - 1 : cds_pos - 1 + len(cds_ref)]
            if cds_pos < 1 or cds_pos + len(cds_ref) - 1 > len(gene_sequence):
                unresolved_rows.append(
                    {
                        "gene": gene,
                        "mutation": mutation,
                        "reason": "c_variant_out_of_bounds",
                    }
                )
                continue
            if observed != cds_ref:
                unresolved_rows.append(
                    {
                        "gene": gene,
                        "mutation": mutation,
                        "reason": f"c_variant_ref_mismatch:{cds_ref}!={observed}",
                    }
                )
                continue
            exact_variants.append((cds_pos, cds_ref, cds_alt))
        else:
            variant_key = f"{gene}_{mutation}"
            mapped_rows = coordinate_rows.get(variant_key, [])
            if not mapped_rows:
                unresolved_rows.append(
                    {
                        "gene": gene,
                        "mutation": mutation,
                        "reason": "missing_coordinate_mapping",
                    }
                )
                continue
            for mapped in mapped_rows:
                try:
                    cds_pos, cds_ref, cds_alt = convert_genomic_variant_to_cds(
                        gene,
                        gene_sequence,
                        gene_features[gene],
                        int(mapped["position"]),
                        str(mapped["ref"]),
                        str(mapped["alt"]),
                    )
                except ValueError as exc:
                    mapping_errors.append(str(exc))
                    continue
                exact_variants.append((cds_pos, cds_ref, cds_alt))
        deduped_variants = sorted(
            set(exact_variants), key=lambda item: (item[0], item[1], item[2])
        )
        if not deduped_variants:
            unresolved_rows.append(
                {
                    "gene": gene,
                    "mutation": mutation,
                    "reason": (
                        mapping_errors[0]
                        if mapping_errors
                        else "no_exact_variants_after_dedup"
                    ),
                }
            )
            continue
        for cds_pos, cds_ref, cds_alt in deduped_variants:
            insert_row = (
                int(gene_record["id"]),
                cds_pos,
                cds_ref,
                cds_alt,
                mutation,
                "CDS",
                None,
                gene_sequence,
                "evo2_7b",
                drug_id,
                "PENDING",
            )
            if insert_row in seen_rows:
                continue
            seen_rows.add(insert_row)
            inserts.append(insert_row)
            realized += 1
            imported_by_gene[gene] += 1
            imported_by_region_type["CDS"] += 1
        if realized == 0:
            unresolved_rows.append(
                {
                    "gene": gene,
                    "mutation": mutation,
                    "reason": "all_exact_variants_deduped_or_invalid",
                }
            )
            continue
        who_row_realizations[f"{gene}:{mutation}"] = realized
    manifest = {
        "generic_rows": generic_rows,
        "unresolved_rows": unresolved_rows,
        "imported_by_gene": dict(imported_by_gene),
        "imported_by_region_type": dict(imported_by_region_type),
        "who_row_realizations": dict(sorted(who_row_realizations.items())),
    }
    return (inserts, manifest)


def main() -> int:
    args = parse_args()
    modeled_genes = resolve_modeled_genes(args.extra_modeled_gene)
    benchmark_root = args.benchmark_root or default_benchmark_root(args.benchmarks_root)
    benchmark_root.mkdir(parents=True, exist_ok=False)
    paths = create_benchmark_layout(benchmark_root)
    shutil.copy2(args.db, paths["db_path"])
    who_rows, gene_window_specs = select_supported_who_rows(
        args.xlsx, args.drug, args.truth, args.config, modeled_genes
    )
    coordinate_rows = load_coordinate_rows(args.xlsx)
    gene_features = load_gene_features(args.gbk, modeled_genes)
    regulatory_sequences = build_regulatory_sequences_by_gene(
        args.config, args.gbk, modeled_genes
    )
    with sqlite3.connect(paths["db_path"]) as conn:
        conn.row_factory = sqlite3.Row
        reset_mutations_table(conn)
        gene_records = load_benchmark_gene_records(conn)
        drug_id = ensure_bedaquiline_id(conn, args.drug)
        insert_rows, manifest = build_exact_variants(
            who_rows,
            coordinate_rows,
            gene_records,
            gene_features,
            regulatory_sequences,
            gene_window_specs,
            drug_id,
        )
        imported_count = import_variants(conn, insert_rows)
        conn.commit()
    write_benchmark_config(args.config, paths["config_path"], paths)
    summary = {
        "drug": args.drug,
        "truth_filter": args.truth,
        "source_db": str(args.db),
        "who_xlsx": str(args.xlsx),
        "genbank": str(args.gbk),
        "benchmark_root": str(paths["root"]),
        "benchmark_db": str(paths["db_path"]),
        "benchmark_config": str(paths["config_path"]),
        "benchmark_export_csv": str(paths["export_csv_path"]),
        "benchmark_analysis_snapshot_csv": str(paths["analysis_snapshot_csv_path"]),
        "who_current_supported_truth_rows": int(len(who_rows)),
        "generic_non_executable_rows": int(len(manifest["generic_rows"])),
        "unresolved_rows": int(len(manifest["unresolved_rows"])),
        "imported_exact_variants": imported_count,
        "imported_by_gene": manifest["imported_by_gene"],
        "imported_by_region_type": manifest["imported_by_region_type"],
        "who_row_realizations": manifest["who_row_realizations"],
        "generic_rows": manifest["generic_rows"],
        "unresolved_details": manifest["unresolved_rows"],
    }
    paths["manifest_path"].write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
