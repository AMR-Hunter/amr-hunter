from __future__ import annotations
import sys
from argparse import ArgumentParser
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple
from urllib.error import URLError
from urllib.request import Request, urlopen
import yaml
from Bio import SeqIO

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
from core.database import DatabaseManager
from core.generator import MutationEngine


def _safe_file_size(path: Path) -> Optional[int]:
    try:
        return path.stat().st_size
    except OSError:
        return None


def _format_size(size_bytes: Optional[int]) -> str:
    if size_bytes is None:
        return "unknown"
    if size_bytes < 1024:
        return f"{size_bytes} B"
    if size_bytes < 1024 * 1024:
        return f"{size_bytes / 1024:.1f} KB"
    return f"{size_bytes / (1024 * 1024):.2f} MB"


def _download_file(url: str, destination: Path, timeout_seconds: int) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    req = Request(url, headers={"User-Agent": "AMR-Hunter/1.0 (+auto-fetch)"})
    with urlopen(req, timeout=max(1, timeout_seconds)) as response:
        content = response.read()
    destination.write_bytes(content)


def _resolve_download_url(
    item_type: str,
    path_str: str,
    species: Optional[str],
    gene: Optional[str],
    fetch_cfg: Dict[str, Any],
) -> Optional[str]:
    path_key = path_str.strip()
    basename = Path(path_key).name
    direct_urls = fetch_cfg.get("direct_urls", {}) or {}
    if path_key in direct_urls:
        return str(direct_urls[path_key])
    if basename in direct_urls:
        return str(direct_urls[basename])
    if item_type == "reference":
        ref_urls = fetch_cfg.get("reference_genomes", {}) or {}
        if species and species in ref_urls:
            return str(ref_urls[species])
    if item_type == "structure":
        structure_urls = fetch_cfg.get("structure_pdb", {}) or {}
        if gene and gene in structure_urls:
            return str(structure_urls[gene])
    if item_type == "ligand":
        ligand_urls = fetch_cfg.get("ligand_sdf", {}) or {}
        if gene and gene in ligand_urls:
            return str(ligand_urls[gene])
        if basename in ligand_urls:
            return str(ligand_urls[basename])
    return None


def _collect_required_assets(config: Dict[str, Any]) -> List[Dict[str, Optional[str]]]:
    assets: List[Dict[str, Optional[str]]] = []
    for ref in config.get("reference_genomes", []) or []:
        path = ref.get("path")
        species = ref.get("species")
        if path:
            assets.append(
                {
                    "item_type": "reference",
                    "path": str(path),
                    "species": str(species) if species else None,
                    "gene": None,
                    "label": f"reference[{species or 'unknown'}]",
                }
            )
    for gene_cfg in config.get("target_genes", []) or []:
        gene = gene_cfg.get("gene")
        structure_pdb = gene_cfg.get("structure_pdb")
        ligand_sdf = gene_cfg.get("ligand_sdf")
        species = gene_cfg.get("species")
        if structure_pdb:
            assets.append(
                {
                    "item_type": "structure",
                    "path": str(structure_pdb),
                    "species": str(species) if species else None,
                    "gene": str(gene) if gene else None,
                    "label": f"structure[{gene or 'unknown'}]",
                }
            )
        if ligand_sdf:
            assets.append(
                {
                    "item_type": "ligand",
                    "path": str(ligand_sdf),
                    "species": str(species) if species else None,
                    "gene": str(gene) if gene else None,
                    "label": f"ligand[{gene or 'unknown'}]",
                }
            )
    unique: Dict[str, Dict[str, Optional[str]]] = {}
    for item in assets:
        path = item.get("path")
        if not path:
            continue
        if path not in unique:
            unique[path] = item
    return list(unique.values())


def ensure_required_assets(config: Dict[str, Any]) -> None:
    fetch_cfg = config.get("data_autofetch", {}) or {}
    enabled = bool(fetch_cfg.get("enabled", False))
    timeout_seconds = int(fetch_cfg.get("timeout_seconds", 60) or 60)
    required_assets = _collect_required_assets(config)
    if not required_assets:
        return
    print("\nChecking required input assets...")
    missing_without_url: List[str] = []
    for asset in required_assets:
        path_str = asset.get("path")
        if not path_str:
            continue
        path = Path(path_str)
        label = asset.get("label") or asset.get("item_type") or "asset"
        if path.exists():
            size_text = _format_size(_safe_file_size(path))
            print(f"  [OK] {label}: {path} ({size_text})")
            continue
        if not enabled:
            missing_without_url.append(f"{label}: {path}")
            continue
        url = _resolve_download_url(
            item_type=str(asset.get("item_type") or ""),
            path_str=path_str,
            species=asset.get("species"),
            gene=asset.get("gene"),
            fetch_cfg=fetch_cfg,
        )
        if not url:
            missing_without_url.append(f"{label}: {path} (no download URL configured)")
            continue
        print(f"  [FETCH] {label}: {url}")
        try:
            _download_file(url, path, timeout_seconds=timeout_seconds)
            size_text = _format_size(_safe_file_size(path))
            print(f"  [DONE]  {label}: {path} ({size_text})")
        except (URLError, OSError) as exc:
            raise RuntimeError(f"Failed downloading {label} from {url}: {exc}") from exc
    if missing_without_url:
        missing_text = "\n".join((f"- {item}" for item in missing_without_url))
        raise FileNotFoundError(
            f"Missing required assets and unable to auto-fetch some files:\n{missing_text}\nPlease provide files manually or configure data_autofetch URLs in config/config.yaml."
        )


def load_config(config_path: Path) -> Dict[str, Any]:
    with config_path.open("r", encoding="utf-8") as file:
        return yaml.safe_load(file) or {}


def upsert_species(conn, name: str, genome_path: Optional[str]) -> int:
    conn.execute(
        "\n        INSERT INTO species (name, genome_path)\n        VALUES (?, ?)\n        ON CONFLICT(name) DO UPDATE SET genome_path=excluded.genome_path\n        ",
        (name, genome_path),
    )
    row = conn.execute("SELECT id FROM species WHERE name = ?", (name,)).fetchone()
    return int(row["id"])


def upsert_drug(conn, name: str, sdf_path: Optional[str] = None) -> int:
    conn.execute(
        "\n        INSERT INTO drugs (name, sdf_path)\n        VALUES (?, ?)\n        ON CONFLICT(name) DO UPDATE SET sdf_path=COALESCE(excluded.sdf_path, drugs.sdf_path)\n        ",
        (name, sdf_path),
    )
    row = conn.execute("SELECT id FROM drugs WHERE name = ?", (name,)).fetchone()
    return int(row["id"])


def upsert_gene(
    conn,
    species_id: int,
    gene_name: str,
    sequence: str,
    logic_type: str = "POSITIVE",
    scenario: Optional[str] = None,
    pathway_group: Optional[str] = None,
    pdb_path: Optional[str] = None,
    dna_sequence: Optional[str] = None,
    ligand_sdf_path: Optional[str] = None,
) -> int:
    conn.execute(
        "\n        INSERT INTO genes \n        (species_id, name, sequence, logic_type, scenario, pathway_group, pdb_path, dna_sequence, ligand_sdf_path)\n        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)\n        ON CONFLICT(species_id, name) DO UPDATE SET\n            sequence=excluded.sequence,\n            logic_type=COALESCE(excluded.logic_type, genes.logic_type),\n            scenario=COALESCE(excluded.scenario, genes.scenario),\n            pathway_group=COALESCE(excluded.pathway_group, genes.pathway_group),\n            pdb_path=COALESCE(excluded.pdb_path, genes.pdb_path),\n            dna_sequence=COALESCE(excluded.dna_sequence, genes.dna_sequence),\n            ligand_sdf_path=COALESCE(excluded.ligand_sdf_path, genes.ligand_sdf_path)\n        ",
        (
            species_id,
            gene_name,
            sequence,
            logic_type,
            scenario,
            pathway_group,
            pdb_path,
            dna_sequence,
            ligand_sdf_path,
        ),
    )
    row = conn.execute(
        "SELECT id FROM genes WHERE species_id = ? AND name = ?",
        (species_id, gene_name),
    ).fetchone()
    return int(row["id"])


def _get_feature_gene_name(feature) -> Optional[str]:
    qualifiers = feature.qualifiers
    if "gene" in qualifiers and qualifiers["gene"]:
        return qualifiers["gene"][0]
    if "locus_tag" in qualifiers and qualifiers["locus_tag"]:
        return qualifiers["locus_tag"][0]
    return None


def _slice_oriented_sequence(record_seq, start: int, end: int, strand: int) -> str:
    if start >= end:
        return ""
    seq = str(record_seq[start:end]).upper().replace("U", "T")
    if not seq:
        return ""
    if strand < 0:
        return str(record_seq[start:end].reverse_complement()).upper().replace("U", "T")
    return seq


def _extract_intergenic_between(
    record_seq,
    anchor_info: Dict[str, Any],
    partner_info: Dict[str, Any],
    trim_left_bp: int = 0,
    trim_right_bp: int = 0,
) -> str:
    anchor_start = int(anchor_info["start"])
    anchor_end = int(anchor_info["end"])
    partner_start = int(partner_info["start"])
    partner_end = int(partner_info["end"])
    if anchor_end <= partner_start:
        inter_start, inter_end = (anchor_end, partner_start)
    elif partner_end <= anchor_start:
        inter_start, inter_end = (partner_end, anchor_start)
    else:
        return ""
    inter_start = inter_start + max(0, trim_left_bp)
    inter_end = inter_end - max(0, trim_right_bp)
    if inter_start >= inter_end:
        return ""
    strand = int(anchor_info.get("strand", 1) or 1)
    return _slice_oriented_sequence(record_seq, inter_start, inter_end, strand)


def _extract_window_sequence(
    record,
    feature_info: Dict[str, Any],
    window: Dict[str, Any],
    gene_features: Dict[str, Dict[str, Any]],
) -> str:
    mode = str(window.get("mode", "upstream")).lower()
    start = int(feature_info["start"])
    end = int(feature_info["end"])
    strand = int(feature_info.get("strand", 1) or 1)
    seq_len = len(record.seq)
    if mode == "upstream":
        length_bp = max(
            1, int(window.get("length_bp", window.get("upstream_bp", 120)) or 120)
        )
        if strand >= 0:
            win_start = max(0, start - length_bp)
            win_end = start
        else:
            win_start = end
            win_end = min(seq_len, end + length_bp)
        return _slice_oriented_sequence(record.seq, win_start, win_end, strand)
    if mode == "downstream":
        length_bp = max(1, int(window.get("length_bp", 60) or 60))
        if strand >= 0:
            win_start = end
            win_end = min(seq_len, end + length_bp)
        else:
            win_start = max(0, start - length_bp)
            win_end = start
        return _slice_oriented_sequence(record.seq, win_start, win_end, strand)
    if mode == "around_start":
        upstream_bp = max(0, int(window.get("upstream_bp", 80) or 80))
        downstream_bp = max(0, int(window.get("downstream_bp", 20) or 20))
        if strand >= 0:
            win_start = max(0, start - upstream_bp)
            win_end = min(seq_len, start + downstream_bp)
        else:
            win_start = max(0, end - downstream_bp)
            win_end = min(seq_len, end + upstream_bp)
        return _slice_oriented_sequence(record.seq, win_start, win_end, strand)
    if mode == "intergenic":
        partner_gene = window.get("partner_gene")
        partner_info = gene_features.get(partner_gene) if partner_gene else None
        if not partner_info:
            return ""
        return _extract_intergenic_between(
            record.seq,
            anchor_info=feature_info,
            partner_info=partner_info,
            trim_left_bp=int(window.get("trim_left_bp", 0) or 0),
            trim_right_bp=int(window.get("trim_right_bp", 0) or 0),
        )
    return ""


def extract_target_gene_contexts(
    genome_path: Path,
    target_genes: Set[str],
    gene_meta_by_name: Dict[str, Dict[str, Any]],
    default_windows: List[Dict[str, Any]],
    intergenic_regions: List[Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    if not genome_path.exists():
        raise FileNotFoundError(f"Reference genome not found: {genome_path}")
    record = SeqIO.read(str(genome_path), "genbank")
    found: Dict[str, Dict[str, Any]] = {}
    gene_features: Dict[str, Dict[str, Any]] = {}
    for feature in record.features:
        if feature.type != "CDS":
            continue
        gene_name = _get_feature_gene_name(feature)
        if not gene_name or gene_name not in target_genes or gene_name in gene_features:
            continue
        location = feature.location
        strand = int(location.strand or 1)
        start = int(location.start)
        end = int(location.end)
        seq = _slice_oriented_sequence(record.seq, start, end, strand)
        if not seq or len(seq) % 3 != 0:
            continue
        gene_features[gene_name] = {
            "start": start,
            "end": end,
            "strand": strand,
            "cds_sequence": seq,
        }
    for gene_name, feature_info in gene_features.items():
        gene_meta = gene_meta_by_name.get(gene_name, {})
        context: Dict[str, Any] = {
            "cds_sequence": feature_info["cds_sequence"],
            "regulatory_regions": [],
        }
        windows: List[Dict[str, Any]] = []
        if gene_meta.get("apply_default_windows", False):
            windows.extend(default_windows)
        windows.extend(gene_meta.get("noncoding_windows", []) or [])
        if gene_meta.get("include_upstream", False):
            windows.append(
                {
                    "name": f"UPSTREAM_{int(gene_meta.get('upstream_bp', 120) or 120)}BP",
                    "mode": "upstream",
                    "length_bp": int(gene_meta.get("upstream_bp", 120) or 120),
                }
            )
        for idx, window in enumerate(windows):
            region_seq = _extract_window_sequence(
                record, feature_info, window, gene_features
            )
            if not region_seq:
                continue
            region_name = window.get("name") or f"REGION_{idx + 1}"
            context["regulatory_regions"].append(
                {"name": str(region_name), "sequence": region_seq}
            )
        found[gene_name] = context
    for region in intergenic_regions:
        anchor_gene = region.get("anchor_gene")
        partner_gene = region.get("partner_gene")
        if not anchor_gene or not partner_gene:
            continue
        if (
            anchor_gene not in found
            or anchor_gene not in gene_features
            or partner_gene not in gene_features
        ):
            continue
        intergenic_seq = _extract_intergenic_between(
            record.seq,
            anchor_info=gene_features[anchor_gene],
            partner_info=gene_features[partner_gene],
            trim_left_bp=int(region.get("trim_left_bp", 0) or 0),
            trim_right_bp=int(region.get("trim_right_bp", 0) or 0),
        )
        if not intergenic_seq:
            continue
        found[anchor_gene]["regulatory_regions"].append(
            {
                "name": str(
                    region.get("name") or f"INTERGENIC_{anchor_gene}_{partner_gene}"
                ),
                "sequence": intergenic_seq,
            }
        )
    return found


def build_task_index(tasks: List[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    index: Dict[str, Dict[str, Any]] = {}
    for task in tasks:
        if not task.get("enabled", True):
            continue
        species = task.get("species")
        gene = task.get("gene")
        if not species or not gene:
            continue
        key = f"{species}::{gene}"
        if key not in index:
            index[key] = {}
        for attr in ["logic_type", "scenario", "pathway_group"]:
            if attr in task:
                index[key][attr] = task[attr]
        drug = task.get("drug")
        if drug:
            if "drugs" not in index[key]:
                index[key]["drugs"] = set()
            index[key]["drugs"].add(drug)
    return index


def main() -> None:
    parser = ArgumentParser(description="Initialize AMR-Hunter database from config")
    parser.add_argument(
        "--config",
        type=Path,
        default=PROJECT_ROOT / "config" / "config.yaml",
        help="Path to config.yaml",
    )
    parser.add_argument(
        "--evo-model",
        default=None,
        help="Evo2 model name to tag generated mutations with (e.g. evo2_7b, evo2_20b). Defaults to the value in config models.evo2.model_name, or 'evo2_7b'.",
    )
    args = parser.parse_args()
    config_path = args.config
    config = load_config(config_path)
    ensure_required_assets(config)
    evo_model: str = args.evo_model or str(
        (config.get("models", {}) or {}).get("evo2", {}).get("model_name") or "evo2_7b"
    )
    local_db_path = config["paths"]["local_db_path"]
    share_db_path = config["paths"]["share_db_path"]
    reference_map = {
        item["species"]: item["path"]
        for item in config.get("reference_genomes", [])
        if item.get("species") and item.get("path")
    }
    target_genes_meta: Dict[str, Dict[str, Any]] = {}
    target_genes_meta_by_species_gene: Dict[str, Dict[str, Any]] = {}
    target_genes_set: Set[str] = set()
    noncoding_cfg = config.get("noncoding", {})
    default_windows = noncoding_cfg.get("default_windows", []) or []
    apply_default_windows = bool(noncoding_cfg.get("apply_default_windows", False))
    intergenic_regions_all = noncoding_cfg.get("intergenic_regions", []) or []
    for item in config.get("target_genes", []):
        species = item.get("species")
        gene = item.get("gene")
        if not species or not gene:
            continue
        key = f"{species}::{gene}"
        target_genes_meta[key] = {
            "logic_type": item.get("logic_type", "POSITIVE"),
            "scenario": item.get("scenario"),
            "pathway_group": item.get("pathway_group"),
            "pdb_path": item.get("structure_pdb"),
            "dna_sequence": item.get("dna_sequence"),
            "ligand_sdf_path": item.get("ligand_sdf"),
            "include_upstream": bool(item.get("include_upstream", False)),
            "upstream_bp": int(item.get("upstream_bp", 0) or 0),
            "noncoding_windows": item.get("noncoding_windows", []) or [],
            "apply_default_windows": bool(
                item.get("apply_default_windows", apply_default_windows)
            ),
        }
        target_genes_meta_by_species_gene[key] = target_genes_meta[key]
        target_genes_set.add(gene)
    task_index = build_task_index(config.get("tasks", []))
    db = DatabaseManager(local_db_path=local_db_path, share_db_path=share_db_path)
    try:
        for species, genome_path in reference_map.items():
            print(f"\nProcessing {species}...")
            species_gene_meta = {
                gene_key.split("::", 1)[1]: meta
                for gene_key, meta in target_genes_meta_by_species_gene.items()
                if gene_key.startswith(f"{species}::")
            }
            species_intergenic = [
                item
                for item in intergenic_regions_all
                if item.get("species") == species
            ]
            contexts = extract_target_gene_contexts(
                Path(genome_path),
                target_genes_set,
                species_gene_meta,
                default_windows,
                species_intergenic,
            )
            print(f"  Found {len(contexts)} target genes in reference")
            species_id = upsert_species(db.conn, species, str(genome_path))
            for gene_name, context in contexts.items():
                sequence = context["cds_sequence"]
                gene_key = f"{species}::{gene_name}"
                gene_meta = target_genes_meta.get(gene_key, {})
                task_meta = task_index.get(gene_key, {})
                logic_type = task_meta.get("logic_type") or gene_meta.get(
                    "logic_type", "POSITIVE"
                )
                scenario = task_meta.get("scenario") or gene_meta.get("scenario")
                pathway_group = task_meta.get("pathway_group") or gene_meta.get(
                    "pathway_group"
                )
                pdb_path = gene_meta.get("pdb_path")
                dna_sequence = gene_meta.get("dna_sequence")
                ligand_sdf_path = gene_meta.get("ligand_sdf_path")
                gene_id = upsert_gene(
                    db.conn,
                    species_id=species_id,
                    gene_name=gene_name,
                    sequence=sequence,
                    logic_type=logic_type,
                    scenario=scenario,
                    pathway_group=pathway_group,
                    pdb_path=pdb_path,
                    dna_sequence=dna_sequence,
                    ligand_sdf_path=ligand_sdf_path,
                )
                print(f"    {gene_name}: {logic_type} - {scenario or 'PROTEIN_LIGAND'}")
                for drug_name in task_meta.get("drugs", set()):
                    drug_id = upsert_drug(db.conn, drug_name)
                engine = MutationEngine(config_path=str(config_path), db_manager=db)
                generated = engine.generate_saturation_mutagenesis(
                    {"gene_id": gene_id, "sequence": sequence, "target_drug_id": None}
                )
                mutation_count = engine.bulk_insert_mutations(
                    generated, evo_model=evo_model
                )
                print(f"      Generated {mutation_count} mutations [{evo_model}]")
                for region in context.get("regulatory_regions", []):
                    regulatory = engine.generate_regulatory_mutagenesis(
                        {
                            "gene_id": gene_id,
                            "sequence": region.get("sequence", ""),
                            "target_drug_id": None,
                        },
                        region_name=str(region.get("name") or "REGULATORY"),
                    )
                    regulatory_count = engine.bulk_insert_mutations(
                        regulatory, evo_model=evo_model
                    )
                    if regulatory_count > 0:
                        print(
                            f"      Generated {regulatory_count} noncoding mutations [{region.get('name')}]"
                        )
            db.conn.commit()
            db.sync_to_share()
        print("\nDatabase initialization complete!")
        print(f"Local DB: {local_db_path}")
        print(f"Shared DB: {share_db_path}")
    finally:
        db.close()


if __name__ == "__main__":
    main()
