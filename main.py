from __future__ import annotations
from collections import deque
import csv
import hashlib
import json
import logging
import os
import subprocess
import sys
import time
from argparse import ArgumentParser
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, as_completed, wait
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple
import yaml
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from core.analyzer import (
    PathwayAnalyzer,
    apply_differential_filtering,
    attach_analysis_decision,
    format_loss_of_function_interpretation,
    has_loss_of_function_signature,
    resolve_analysis_decision,
)
from core.clients import Boltz2Client, Evo2Client
from core.database import DatabaseManager
from core.generator import MutationEngine

FAST_RUNTIME_BOLTZ_ARTIFACTS_ROOT = Path(
    "/root/gpufree-data/amr_hunter_runtime/results/boltz_artifacts"
)
SHARE_ROOT = Path("/root/gpufree-share")
ANALYSIS_EXPORT_HEADER = [
    "gene",
    "logic_type",
    "scenario",
    "mutation",
    "evo_model",
    "evo_delta",
    "stability_score",
    "ptm",
    "iptm",
    "plddt",
    "boltz_metric_context",
    "aux_binding_affinity",
    "aux_affinity_pred_value",
    "aux_pair_energy",
    "aux_complex_energy",
    "synergy_score",
    "synergy_tags",
    "resistance_label",
    "functional_state",
    "evidence_code",
    "final_status",
    "interpretation",
]


def load_config(config_path: Path) -> Dict[str, Any]:
    with config_path.open("r", encoding="utf-8") as file:
        return yaml.safe_load(file) or {}


def normalize_gene_filter(raw_genes: Optional[List[str]]) -> Optional[List[str]]:
    if not raw_genes:
        return None
    genes: List[str] = []
    seen: Set[str] = set()
    for raw_gene in raw_genes:
        gene = str(raw_gene).strip()
        if not gene or gene in seen:
            continue
        genes.append(gene)
        seen.add(gene)
    return genes or None


def load_mutation_id_filter(path: Optional[Path]) -> Optional[Set[int]]:
    if path is None:
        return None
    mutation_ids: Set[int] = set()
    with path.open("r", encoding="utf-8") as handle:
        for line_no, raw_line in enumerate(handle, start=1):
            line = raw_line.split("#", 1)[0].strip()
            if not line:
                continue
            for token in line.replace(",", " ").split():
                try:
                    mutation_id = int(token)
                except ValueError as exc:
                    raise ValueError(
                        f"{path}:{line_no}: invalid mutation id {token!r}"
                    ) from exc
                if mutation_id <= 0:
                    raise ValueError(
                        f"{path}:{line_no}: mutation id must be positive: {mutation_id}"
                    )
                mutation_ids.add(mutation_id)
    if not mutation_ids:
        raise ValueError(f"{path}: mutation id file did not contain any ids")
    return mutation_ids


def build_target_filter_sql(
    gene_names: Optional[List[str]], mutation_ids: Optional[Set[int]]
) -> Tuple[str, List[Any]]:
    clauses: List[str] = []
    params: List[Any] = []
    if gene_names:
        placeholders = ", ".join(("?" for _ in gene_names))
        clauses.append(f"g.name IN ({placeholders})")
        params.extend(gene_names)
    if mutation_ids:
        ordered_ids = sorted(mutation_ids)
        placeholders = ", ".join(("?" for _ in ordered_ids))
        clauses.append(f"m.id IN ({placeholders})")
        params.extend(ordered_ids)
    if not clauses:
        return ("", [])
    return (" AND " + " AND ".join((f"({clause})" for clause in clauses)), params)


def setup_logger(log_path: Path) -> logging.Logger:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("amr_hunter.pipeline")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    file_handler = logging.FileHandler(log_path, encoding="utf-8")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)
    return logger


def append_batch_metrics(metrics_path: Path, payload: Dict[str, Any]) -> None:
    metrics_path.parent.mkdir(parents=True, exist_ok=True)
    event = {"timestamp": datetime.now(timezone.utc).isoformat(), **payload}
    with metrics_path.open("a", encoding="utf-8") as fp:
        fp.write(json.dumps(event, ensure_ascii=False) + "\n")


def _default_boltz_artifacts_root(export_path: Path) -> Path:
    candidate = export_path.parent / "boltz_artifacts"
    try:
        export_parent = export_path.parent.resolve()
        share_root = SHARE_ROOT.resolve()
        if export_parent == share_root or share_root in export_parent.parents:
            return FAST_RUNTIME_BOLTZ_ARTIFACTS_ROOT
    except OSError:
        pass
    return candidate


def _format_analysis_export_row(row: Any) -> Tuple[Tuple[str, ...], List[Any]]:
    resolved = dict(row)
    decision = resolve_analysis_decision(resolved)
    final_status = (
        "SENSITIVE" if resolved["status"] == "RE_SENSITIZED" else resolved["status"]
    )
    scenario = str(resolved.get("scenario") or "").strip().upper()
    has_boltz_metrics = any(
        (
            resolved.get(field) is not None
            for field in (
                "boltz_ptm",
                "boltz_iptm",
                "boltz_plddt",
                "boltz_stability_score",
            )
        )
    )
    folding_only_metrics = scenario in {"FOLDING_ONLY", "PROTEIN_ONLY"}
    if has_boltz_metrics and folding_only_metrics:
        boltz_metric_context = "folding_only_plddt_primary"
        export_ptm = ""
        export_iptm = ""
    elif has_boltz_metrics:
        boltz_metric_context = "ptm_iptm_plddt"
        export_ptm = resolved["boltz_ptm"] if resolved["boltz_ptm"] is not None else ""
        export_iptm = (
            resolved["boltz_iptm"] if resolved["boltz_iptm"] is not None else ""
        )
    else:
        boltz_metric_context = ""
        export_ptm = ""
        export_iptm = ""
    mutation_label = resolved["aa_change"]
    if not mutation_label:
        region_name = resolved["region_name"] or resolved["region_type"] or "NONCODING"
        mutation_label = (
            f"{region_name}:{resolved['ref']}{resolved['pos']}{resolved['alt']}"
        )
    export_row = [
        resolved["name"],
        resolved["logic_type"],
        scenario,
        mutation_label,
        resolved["evo_model"] or "evo2_7b",
        resolved["evo_delta"] or "",
        resolved["boltz_stability_score"] or "",
        export_ptm,
        export_iptm,
        resolved["boltz_plddt"] or "",
        boltz_metric_context,
        resolved["boltz_binding_affinity"] or "",
        resolved["boltz_affinity_pred_value"] or "",
        resolved["boltz_pair_energy"] or "",
        resolved["boltz_complex_energy"] or "",
        resolved["synergy_score"] or 0,
        resolved["synergy_tags"] or "",
        decision.resistance_phenotype,
        decision.functional_state,
        decision.evidence_code,
        final_status,
        resolved["final_interpretation"] or "",
    ]
    sort_key = tuple((str(value) for value in export_row))
    return (sort_key, export_row)


def _write_analysis_csv(
    path: Path, results: List[Any], logger: logging.Logger, label: str
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    formatted_rows = [_format_analysis_export_row(row) for row in results]
    formatted_rows.sort(key=lambda item: item[0])
    with path.open("w", encoding="utf-8", newline="") as file_obj:
        writer = csv.writer(file_obj)
        writer.writerow(ANALYSIS_EXPORT_HEADER)
        for _, export_row in formatted_rows:
            writer.writerow(export_row)
    logger.info("%s exported to %s", label, path)


def maybe_checkpoint_sync(
    db: DatabaseManager,
    logger: logging.Logger,
    stage: str,
    processed: int,
    failed: int,
    last_sync_monotonic: float,
    interval_seconds: int,
    *,
    force: bool = False,
) -> float:
    now = time.monotonic()
    if not force and now - last_sync_monotonic < interval_seconds:
        return last_sync_monotonic
    db.conn.commit()
    db.sync_to_share_async()
    logger.info(
        "[检查点] 已触发数据库异步同步到共享盘: stage=%s processed=%d failed=%d",
        stage,
        processed,
        failed,
    )
    return now


def _format_metric_summary(
    metrics: Dict[str, Any], field_names: Tuple[str, ...]
) -> str:
    parts: List[str] = []
    for field_name in field_names:
        value = metrics.get(field_name)
        if value is None:
            continue
        parts.append(f"{field_name}={value}")
    return ", ".join(parts) if parts else "none"


def _resolve_runtime_boltz_context(
    gene_settings: Dict[str, Dict[str, Any]], mut: Dict[str, Any]
) -> Tuple[Dict[str, Any], str, str]:
    gene_name = str(mut.get("gene_name") or mut.get("name") or "")
    runtime = gene_settings.get(gene_name, {})
    scenario = str(
        runtime.get("scenario") or mut.get("scenario") or "PROTEIN_LIGAND"
    ).upper()
    region_type = str(mut.get("region_type") or "CDS").upper()
    return (runtime, scenario, region_type)


def _score_boltz_candidate(
    config: Dict[str, Any],
    gene_settings: Dict[str, Dict[str, Any]],
    mut: Dict[str, Any],
    evo_delta: float,
) -> Tuple[Optional[float], str, str]:
    runtime, scenario, region_type = _resolve_runtime_boltz_context(gene_settings, mut)
    logic_type = str(mut.get("logic_type") or "POSITIVE").upper()
    evo_cfg = (config.get("pipeline") or {}).get("evo") or {}
    lethal_threshold = float(evo_cfg.get("lethal_threshold", -10.0))
    collapse_threshold = float(evo_cfg.get("collapse_threshold", -8.0))
    if region_type == "REGULATORY" and scenario != "PROTEIN_DNA":
        return (
            None,
            "regulatory candidate without a protein-DNA Boltz scenario",
            "ineligible",
        )
    if bool(runtime.get("force_boltz", False)):
        score = 700.0 - abs(evo_delta + 3.0) * 15.0
        return (
            score,
            "force_boltz target reserved for high-value exploration",
            "force_boltz",
        )
    if logic_type == "NEGATIVE":
        target = collapse_threshold + 0.5
        score = 520.0 - abs(evo_delta - target) * 42.0
        if evo_delta >= 0:
            score -= 120.0 + evo_delta * 8.0
        bucket = "resistance"
        reason = f"NEGATIVE candidate prioritized by closeness to LOF threshold ({collapse_threshold:.1f})"
    elif logic_type == "STRUCTURAL":
        target = collapse_threshold + 1.5
        score = 360.0 - abs(evo_delta - target) * 32.0
        if evo_delta < collapse_threshold:
            score -= 60.0
        bucket = "resistance"
        reason = f"STRUCTURAL candidate prioritized near partial-collapse band above {collapse_threshold:.1f}"
    else:
        target = -3.0
        score = 220.0 - abs(evo_delta - target) * 30.0
        if evo_delta < lethal_threshold:
            score -= 120.0
        bucket = "exploration"
        reason = f"POSITIVE candidate kept in exploration band around {target:.1f}"
    if region_type == "REGULATORY":
        score += 10.0
    if scenario == "PROTEIN_DNA":
        score += 5.0
    return (score, reason, bucket)


def _interleave_boltz_candidates(
    force_candidates: List[Tuple[float, Any, str, str]],
    resistance_candidates: List[Tuple[float, Any, str, str]],
    exploration_candidates: List[Tuple[float, Any, str, str]],
    limit: int,
) -> List[Tuple[float, Any, str, str]]:
    selected: List[Tuple[float, Any, str, str]] = []
    queues = {
        "resistance": list(resistance_candidates),
        "force_boltz": list(force_candidates),
        "exploration": list(exploration_candidates),
    }
    cycle = ["resistance", "resistance", "resistance", "force_boltz", "exploration"]
    while len(selected) < limit and any(queues.values()):
        progressed = False
        for bucket in cycle:
            if len(selected) >= limit:
                break
            queue = queues[bucket]
            if not queue:
                continue
            selected.append(queue.pop(0))
            progressed = True
        if progressed:
            continue
        for fallback_bucket in ("resistance", "force_boltz", "exploration"):
            if len(selected) >= limit:
                break
            queue = queues[fallback_bucket]
            if not queue:
                continue
            selected.append(queue.pop(0))
    return selected


def build_runtime_gene_settings(config: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    settings: Dict[str, Dict[str, Any]] = {}
    for section_name in ("target_genes", "tasks"):
        for item in config.get(section_name, []) or []:
            if section_name == "tasks" and (not item.get("enabled", True)):
                continue
            gene_name = item.get("gene")
            if not gene_name:
                continue
            gene_settings = settings.setdefault(str(gene_name), {})
            for key in ("scenario", "dna_sequence", "complex_partners", "force_boltz"):
                value = item.get(key)
                if value is None:
                    continue
                if isinstance(value, list) and (not value):
                    continue
                gene_settings[key] = value
    return settings


def apply_single_base_mutation(
    sequence: str, pos: int, ref: Optional[str], alt: str
) -> str:
    normalized = (sequence or "").upper().replace("U", "T")
    if not normalized:
        raise ValueError("Missing source sequence for mutation application")
    if pos < 1 or pos > len(normalized):
        raise ValueError(
            f"Mutation position {pos} out of range for sequence length {len(normalized)}"
        )
    expected = (ref or "").upper().replace("U", "T")
    replacement = (alt or "").upper().replace("U", "T")
    start_index = pos - 1
    end_index = start_index + len(expected)
    if not expected and (not replacement):
        raise ValueError("Mutation edit must provide a non-empty ref or alt")
    if end_index > len(normalized):
        raise ValueError(
            f"Mutation span [{pos}, {pos + len(expected) - 1}] exceeds sequence length {len(normalized)}"
        )
    observed = normalized[start_index:end_index] if expected else ""
    if expected and observed != expected:
        raise ValueError(
            f"Reference base mismatch at position {pos}: expected {expected}, observed {observed}"
        )
    return normalized[:start_index] + replacement + normalized[end_index:]


def translate_coding_sequence(sequence: str) -> str:
    normalized = (sequence or "").upper().replace("U", "T")
    if not normalized:
        raise ValueError("Missing coding sequence for translation")
    if len(normalized) % 3 != 0:
        raise ValueError(
            f"Coding sequence length {len(normalized)} is not divisible by 3; cannot translate"
        )
    translated = str(Seq(normalized).translate(table=11, to_stop=False))
    if "*" in translated:
        translated = translated.split("*", 1)[0]
    if not translated:
        raise ValueError("Translated protein sequence is empty")
    return translated


def load_partner_sequences(
    db: DatabaseManager,
    partner_names: List[str],
    cache: Dict[Tuple[str, ...], List[str]],
) -> List[str]:
    cache_key = tuple(partner_names)
    if cache_key in cache:
        return cache[cache_key]
    if not partner_names:
        cache[cache_key] = []
        return []
    placeholders = ",".join(("?" for _ in partner_names))
    rows = db.conn.execute(
        f"SELECT name, sequence FROM genes WHERE name IN ({placeholders})",
        tuple(partner_names),
    ).fetchall()
    sequence_by_name = {str(row["name"]): str(row["sequence"] or "") for row in rows}
    missing = [name for name in partner_names if not sequence_by_name.get(name)]
    if missing:
        raise ValueError(f"Missing partner sequences for genes: {', '.join(missing)}")
    sequences = [
        translate_coding_sequence(sequence_by_name[name]) for name in partner_names
    ]
    cache[cache_key] = sequences
    return sequences


def run_evo_stage(
    config: Dict[str, Any],
    db: DatabaseManager,
    logger: logging.Logger,
    evo_model: str = "evo2_7b",
    limit: Optional[int] = None,
    gene_names: Optional[List[str]] = None,
    mutation_ids: Optional[Set[int]] = None,
) -> int:
    logger.info("=" * 80)
    logger.info("STAGE 1: EVO2 FITNESS SCORING  [model=%s]", evo_model)
    logger.info("=" * 80)
    evo_cfg = config["pipeline"]["evo"]
    gene_settings = build_runtime_gene_settings(config)
    evo_client = Evo2Client(
        base_url=config["models"]["evo2"]["base_url"],
        timeout_seconds=config["models"]["evo2"]["timeout_seconds"],
        max_retries=config["models"]["evo2"]["max_retries"],
    )
    limit_clause = f"LIMIT {int(limit)}" if limit else ""
    target_filter_sql, target_filter_params = build_target_filter_sql(
        gene_names, mutation_ids
    )
    if gene_names or mutation_ids:
        logger.info(
            "Target filter active for Evo: genes=%s mutation_ids=%s",
            ",".join(gene_names or ["*"]),
            len(mutation_ids or set()),
        )
    pending = db.conn.execute(
        f"\n        SELECT m.id, m.gene_id, m.pos, m.ref, m.alt, m.aa_change,\n               m.region_type, m.region_name, m.source_sequence,\n             g.name, g.sequence, g.logic_type, g.scenario\n        FROM mutations m\n        JOIN genes g ON m.gene_id = g.id\n        WHERE m.status = 'PENDING' AND m.evo_model = ?\n        {target_filter_sql}\n        ORDER BY g.id, m.pos\n        {limit_clause}\n        ",
        (evo_model, *target_filter_params),
    ).fetchall()
    logger.info(f"Found {len(pending)} PENDING mutations")
    if not pending:
        logger.info("No mutations to process. Stage 1 complete.")
        return 0

    def resolve_evo_status(mut: Dict[str, Any], evo_delta: float) -> Dict[str, Any]:
        candidate = {
            "id": mut["id"],
            "name": mut["gene_name"],
            "gene_name": mut["gene_name"],
            "logic_type": mut["logic_type"],
            "region_type": mut["region_type"],
            "scenario": mut.get("scenario"),
            "evo_delta": evo_delta,
            "aa_change": mut.get("aa_change"),
            "ref": mut.get("ref"),
            "alt": mut.get("alt"),
        }
        resolved = apply_differential_filtering([candidate], evo_cfg, gene_settings)[0]
        _, scenario, region_type = _resolve_runtime_boltz_context(
            gene_settings, resolved
        )
        if (
            region_type == "REGULATORY"
            and resolved.get("status") == "BOLTZ_READY"
            and (scenario != "PROTEIN_DNA")
        ):
            resolved["status"] = "COMPLETED"
            resolved["final_interpretation"] = (
                f"Evo-retained regulatory candidate (delta={evo_delta:.2f}); no protein-DNA Boltz scenario available, kept as Evo-only result"
            )
            attach_analysis_decision(resolved)
            return resolved
        if (
            region_type == "REGULATORY"
            and resolved.get("status") == "BOLTZ_READY"
            and (scenario == "PROTEIN_DNA")
        ):
            resolved["status"] = "BOLTZ_READY"
            resolved["final_interpretation"] = (
                f"Regulatory protein-DNA candidate retained (delta={evo_delta:.2f}); queued for Boltz priority ranking"
            )
            attach_analysis_decision(resolved)
            return resolved
        if resolved.get("status") == "BOLTZ_READY" and region_type == "CDS":
            resolved["status"] = "BOLTZ_READY"
            resolved["final_interpretation"] = (
                f"Evo-retained CDS candidate (delta={evo_delta:.2f}); queued for Boltz priority ranking"
            )
            attach_analysis_decision(resolved)
            return resolved
        attach_analysis_decision(resolved)
        return resolved

    wt_cache: Dict[Tuple[int, str, str], float] = {}
    mutations_by_context: Dict[Tuple[int, str, str], List[Dict[str, Any]]] = {}
    for row in pending:
        gene_id = row["gene_id"]
        region_type = row["region_type"] or "CDS"
        region_name = row["region_name"] or ""
        wt_region_seq = row["source_sequence"] or row["sequence"] or ""
        context_key = (gene_id, region_type, region_name)
        if context_key not in mutations_by_context:
            mutations_by_context[context_key] = []
        mutations_by_context[context_key].append(
            {
                "id": row["id"],
                "gene_id": gene_id,
                "gene_name": row["name"],
                "logic_type": row["logic_type"],
                "sequence": row["sequence"],
                "region_type": row["region_type"] or "CDS",
                "region_name": row["region_name"],
                "source_sequence": wt_region_seq,
                "scenario": gene_settings.get(str(row["name"]), {}).get("scenario")
                or row["scenario"],
                "pos": row["pos"],
                "ref": row["ref"],
                "alt": row["alt"],
                "aa_change": row["aa_change"],
            }
        )
    processed = 0
    failed = 0
    workers = evo_cfg.get("workers", 8)
    checkpoint_interval_seconds = int(
        config.get("pipeline", {}).get("checkpoint_sync_interval_seconds", 300)
    )
    last_checkpoint_sync = time.monotonic()
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {}
        for context_key, muts in mutations_by_context.items():
            if context_key not in wt_cache:
                try:
                    wt_seq = muts[0].get("source_sequence") or muts[0].get("sequence")
                    wt_cache[context_key] = evo_client.get_likelihood(wt_seq)
                    logger.debug(
                        f"Gene {muts[0]['gene_name']} WT score: {wt_cache[context_key]:.4f}"
                    )
                except Exception as exc:
                    logger.warning(
                        f"Failed to score WT for context {context_key}: {exc}"
                    )
                    wt_cache[context_key] = 0.0
            wt_score = wt_cache[context_key]
            for mut in muts:
                if mut["sequence"] and mut["pos"] and (mut.get("alt") is not None):
                    try:
                        wt_region_seq = mut.get("source_sequence") or mut.get(
                            "sequence"
                        )
                        if not wt_region_seq:
                            raise ValueError(
                                "Missing source sequence for mutation scoring"
                            )
                        mt_seq = apply_single_base_mutation(
                            str(wt_region_seq),
                            int(mut["pos"]),
                            mut.get("ref"),
                            str(mut["alt"] or ""),
                        )
                        future = executor.submit(evo_client.get_likelihood, mt_seq)
                        futures[future] = (mut, wt_score)
                    except Exception as exc:
                        logger.error(f"Mutation {mut['id']} prep failed: {exc}")
                        db.conn.execute(
                            "UPDATE mutations SET last_error = ?, retry_count = retry_count + 1, updated_at = CURRENT_TIMESTAMP WHERE id = ?",
                            (str(exc), mut["id"]),
                        )
                        failed += 1
        for future in as_completed(futures):
            mut, wt_score = futures[future]
            mut_id = mut["id"]
            try:
                mt_score = future.result()
                evo_delta = mt_score - wt_score
                resolved = resolve_evo_status(mut, evo_delta)
                status = str(resolved.get("status") or "PENDING")
                final_interpretation = resolved.get("final_interpretation")
                db.conn.execute(
                    "\n                    UPDATE mutations\n                    SET evo_mt_score = ?, evo_delta = ?, status = ?, final_interpretation = ?,\n                        resistance_phenotype = ?, functional_state = ?, evidence_code = ?,\n                        last_error = NULL, updated_at = CURRENT_TIMESTAMP\n                    WHERE id = ?\n                    ",
                    (
                        mt_score,
                        evo_delta,
                        status,
                        final_interpretation,
                        resolved.get("resistance_phenotype"),
                        resolved.get("functional_state"),
                        resolved.get("evidence_code"),
                        mut_id,
                    ),
                )
                processed += 1
                last_checkpoint_sync = maybe_checkpoint_sync(
                    db,
                    logger,
                    "evo",
                    processed,
                    failed,
                    last_checkpoint_sync,
                    checkpoint_interval_seconds,
                )
                if processed % 50 == 0:
                    logger.info(f"Processed {processed} mutations...")
            except Exception as exc:
                logger.error(f"Mutation {mut_id} scoring failed: {exc}")
                db.conn.execute(
                    "UPDATE mutations SET last_error = ?, retry_count = retry_count + 1, updated_at = CURRENT_TIMESTAMP WHERE id = ?",
                    (str(exc), mut_id),
                )
                failed += 1
                last_checkpoint_sync = maybe_checkpoint_sync(
                    db,
                    logger,
                    "evo",
                    processed,
                    failed,
                    last_checkpoint_sync,
                    checkpoint_interval_seconds,
                )
    last_checkpoint_sync = maybe_checkpoint_sync(
        db,
        logger,
        "evo",
        processed,
        failed,
        last_checkpoint_sync,
        checkpoint_interval_seconds,
        force=True,
    )
    logger.info(f"Stage 1 complete: {processed} processed, {failed} failed")
    append_batch_metrics(
        Path(config["pipeline"]["metrics_path"]),
        {"stage": "evo", "processed": processed, "failed": failed},
    )
    return 0 if failed == 0 else 1


def _boltz_cache_key(
    protein_sequence: str, scenario: str, ligand_sdf: Optional[str]
) -> Tuple[str, str, str]:
    seq_hash = hashlib.md5(protein_sequence.encode()).hexdigest()
    ligand_content = ""
    if ligand_sdf:
        try:
            ligand_content = Path(ligand_sdf).read_text(
                encoding="utf-8", errors="replace"
            )
        except OSError:
            ligand_content = ligand_sdf
    ligand_hash = (
        hashlib.md5(ligand_content.encode()).hexdigest() if ligand_content else ""
    )
    return (seq_hash, scenario.upper(), ligand_hash)


def _boltz_cache_lookup(
    db: DatabaseManager, seq_hash: str, scenario: str, ligand_hash: str
) -> Optional[Dict[str, Optional[float]]]:
    row = db.conn.execute(
        "\n        SELECT iptm, ptm, plddt, complex_energy, binding_affinity, affinity_pred_value,\n               pair_energy, stability_score, clash_score\n        FROM boltz_cache\n        WHERE seq_hash = ? AND scenario = ? AND ligand_hash = ?\n        ",
        (seq_hash, scenario, ligand_hash),
    ).fetchone()
    if row is None:
        return None
    metrics = {
        "iptm": row["iptm"],
        "ptm": row["ptm"],
        "plddt": row["plddt"],
        "complex_energy": row["complex_energy"],
        "binding_affinity": row["binding_affinity"],
        "affinity_pred_value": row["affinity_pred_value"],
        "pair_energy": row["pair_energy"],
        "stability_score": row["stability_score"],
        "clash_score": row["clash_score"],
    }
    return None if _boltz_metrics_empty(metrics) else metrics


def _boltz_metrics_empty(metrics: Dict[str, Optional[float]]) -> bool:
    return all(
        (
            metrics.get(key) is None
            for key in (
                "iptm",
                "ptm",
                "plddt",
                "complex_energy",
                "binding_affinity",
                "affinity_pred_value",
                "pair_energy",
                "stability_score",
                "clash_score",
            )
        )
    )


def _boltz_cache_store(
    db: DatabaseManager,
    seq_hash: str,
    scenario: str,
    ligand_hash: str,
    metrics: Dict[str, Optional[float]],
) -> None:
    if _boltz_metrics_empty(metrics):
        return
    db.conn.execute(
        "\n        INSERT OR REPLACE INTO boltz_cache\n        (seq_hash, scenario, ligand_hash, iptm, ptm, plddt, complex_energy, binding_affinity,\n         affinity_pred_value, pair_energy, stability_score, clash_score)\n        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)\n        ",
        (
            seq_hash,
            scenario,
            ligand_hash,
            metrics.get("iptm"),
            metrics.get("ptm"),
            metrics.get("plddt"),
            metrics.get("complex_energy"),
            metrics.get("binding_affinity"),
            metrics.get("affinity_pred_value"),
            metrics.get("pair_energy"),
            metrics.get("stability_score"),
            metrics.get("clash_score"),
        ),
    )


def _boltz_cache_prefetch(
    db: DatabaseManager, cache_keys: List[Tuple[str, str, str]]
) -> Dict[Tuple[str, str, str], Dict[str, Optional[float]]]:
    if not cache_keys:
        return {}
    unique_keys = list(dict.fromkeys(cache_keys))
    clauses = [
        "(seq_hash = ? AND scenario = ? AND ligand_hash = ?)" for _ in unique_keys
    ]
    params: List[str] = []
    for seq_hash, scenario, ligand_hash in unique_keys:
        params.extend((seq_hash, scenario, ligand_hash))
    rows = db.conn.execute(
        f"\n        SELECT seq_hash, scenario, ligand_hash, iptm, ptm, plddt, complex_energy,\n               binding_affinity, affinity_pred_value, pair_energy, stability_score, clash_score\n        FROM boltz_cache\n        WHERE {' OR '.join(clauses)}\n        ",
        params,
    ).fetchall()
    cached: Dict[Tuple[str, str, str], Dict[str, Optional[float]]] = {}
    for row in rows:
        metrics = {
            "iptm": row["iptm"],
            "ptm": row["ptm"],
            "plddt": row["plddt"],
            "complex_energy": row["complex_energy"],
            "binding_affinity": row["binding_affinity"],
            "affinity_pred_value": row["affinity_pred_value"],
            "pair_energy": row["pair_energy"],
            "stability_score": row["stability_score"],
            "clash_score": row["clash_score"],
        }
        if _boltz_metrics_empty(metrics):
            continue
        cached[row["seq_hash"], row["scenario"], row["ligand_hash"]] = metrics
    return cached


def _write_boltz_metrics_to_mutation(
    db: DatabaseManager, mutation_id: int, metrics: Dict[str, Optional[float]]
) -> None:
    db.conn.execute(
        "\n        UPDATE mutations\n        SET status = 'COMPLETED',\n            boltz_iptm = ?,\n            boltz_ptm = ?,\n            boltz_plddt = ?,\n            boltz_complex_energy = ?,\n            boltz_binding_affinity = ?,\n            boltz_affinity_pred_value = ?,\n            boltz_pair_energy = ?,\n            boltz_stability_score = ?,\n            boltz_clash = ?,\n            boltz_artifact_dir = ?,\n            resistance_phenotype = NULL,\n            functional_state = NULL,\n            evidence_code = NULL,\n            last_error = NULL,\n            updated_at = CURRENT_TIMESTAMP\n        WHERE id = ?\n        ",
        (
            metrics.get("iptm"),
            metrics.get("ptm"),
            metrics.get("plddt"),
            metrics.get("complex_energy"),
            metrics.get("binding_affinity"),
            metrics.get("affinity_pred_value"),
            metrics.get("pair_energy"),
            metrics.get("stability_score"),
            metrics.get("clash_score"),
            metrics.get("artifact_dir"),
            mutation_id,
        ),
    )


def _query_gpu_free_memory_gb() -> Optional[float]:
    try:
        result = subprocess.run(
            ["nvidia-smi", "--query-gpu=memory.free", "--format=csv,noheader,nounits"],
            capture_output=True,
            text=True,
            timeout=3,
            check=True,
        )
    except (FileNotFoundError, OSError, subprocess.SubprocessError):
        return None
    values: List[float] = []
    for line in result.stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            values.append(float(line) / 1024.0)
        except ValueError:
            continue
    if not values:
        return None
    return max(values)


def _estimate_boltz_job_vram_gb(
    job: Dict[str, Any], dynamic_cfg: Dict[str, Any]
) -> float:
    defaults = {
        "FOLDING_ONLY": 10.0,
        "PROTEIN_ONLY": 12.0,
        "PROTEIN_LIGAND": 14.0,
        "PROTEIN_DNA": 16.0,
        "DIMER_DNA": 20.0,
        "PROTEIN_PROTEIN": 18.0,
        "COMPLEX": 18.0,
        "DIMER": 16.0,
    }
    configured = dynamic_cfg.get("scenario_vram_gb") or {}
    scenario = str(job["scenario"] or "PROTEIN_LIGAND").upper()
    estimated = float(configured.get(scenario, defaults.get(scenario, 28.0)))
    total_length = len(str(job["protein_sequence"] or ""))
    total_length += sum(
        (len(str(seq or "")) for seq in job.get("partner_sequences") or [])
    )
    total_length += len(str(job.get("dna_sequence") or ""))
    if total_length >= 1200:
        estimated += 5.0
    elif total_length >= 900:
        estimated += 3.0
    elif total_length >= 600:
        estimated += 1.5
    return estimated


def _current_boltz_memory_budget_gb(dynamic_cfg: Dict[str, Any]) -> float:
    configured_budget = float(dynamic_cfg.get("memory_budget_gb", 40.0))
    reserve_gb = float(dynamic_cfg.get("reserve_gb", 6.0))
    free_gb = _query_gpu_free_memory_gb()
    if free_gb is None:
        return configured_budget
    return max(0.0, min(configured_budget, free_gb - reserve_gb))


def _query_gpu_runtime_snapshot() -> Optional[Dict[str, float]]:
    try:
        result = subprocess.run(
            [
                "nvidia-smi",
                "--query-gpu=utilization.gpu,memory.used,memory.total,memory.free",
                "--format=csv,noheader,nounits",
            ],
            capture_output=True,
            text=True,
            timeout=3,
            check=True,
        )
    except (FileNotFoundError, OSError, subprocess.SubprocessError):
        return None
    best: Optional[Dict[str, float]] = None
    for line in result.stdout.splitlines():
        parts = [part.strip() for part in line.split(",")]
        if len(parts) != 4:
            continue
        try:
            utilization = float(parts[0])
            memory_used_gb = float(parts[1]) / 1024.0
            memory_total_gb = float(parts[2]) / 1024.0
            memory_free_gb = float(parts[3]) / 1024.0
        except ValueError:
            continue
        snapshot = {
            "utilization_gpu": utilization,
            "memory_used_gb": memory_used_gb,
            "memory_total_gb": memory_total_gb,
            "memory_free_gb": memory_free_gb,
            "memory_utilization": (
                memory_used_gb / memory_total_gb if memory_total_gb > 0 else 0.0
            ),
        }
        if best is None or snapshot["utilization_gpu"] > best["utilization_gpu"]:
            best = snapshot
    return best


def _boltz_job_total_sequence_length(job: Dict[str, Any]) -> int:
    total_length = len(str(job.get("protein_sequence") or ""))
    total_length += sum(
        (len(str(seq or "")) for seq in job.get("partner_sequences") or [])
    )
    total_length += len(str(job.get("dna_sequence") or ""))
    return total_length


def _estimate_boltz_job_duration_score(job: Dict[str, Any]) -> float:
    scenario = str(job.get("scenario") or "PROTEIN_LIGAND").upper()
    scenario_bias = {
        "FOLDING_ONLY": 0.8,
        "PROTEIN_ONLY": 1.0,
        "PROTEIN_LIGAND": 1.3,
        "PROTEIN_PROTEIN": 1.6,
        "COMPLEX": 1.8,
        "DIMER": 1.5,
        "PROTEIN_DNA": 1.4,
    }
    total_length = _boltz_job_total_sequence_length(job)
    estimated_vram_gb = float(job.get("estimated_vram_gb") or 0.0)
    partner_count = len(job.get("partner_sequences") or [])
    ligand_penalty = 0.5 if job.get("ligand_sdf") else 0.0
    dna_penalty = 0.4 if job.get("dna_sequence") else 0.0
    return (
        total_length / 100.0
        + estimated_vram_gb * 1.6
        + partner_count * 2.0
        + ligand_penalty
        + dna_penalty
        + scenario_bias.get(scenario, 1.1)
    )


def _rebalance_pending_boltz_jobs(
    pending_jobs: List[Dict[str, Any]], dynamic_cfg: Dict[str, Any]
) -> List[Dict[str, Any]]:
    if not pending_jobs or not dynamic_cfg.get("duration_pairing_enabled", True):
        return pending_jobs
    if int(dynamic_cfg.get("max_workers", 1)) <= 1:
        return pending_jobs
    parallel_safe: List[Tuple[float, int, Dict[str, Any]]] = []
    parallel_unsafe: List[Tuple[float, int, Dict[str, Any]]] = []
    for index, job in enumerate(pending_jobs):
        duration_score = _estimate_boltz_job_duration_score(job)
        job["estimated_duration_score"] = duration_score
        item = (duration_score, index, job)
        if _adaptive_parallel_job_eligible(job, dynamic_cfg):
            parallel_safe.append(item)
        else:
            parallel_unsafe.append(item)
    used_keys = set()
    reordered: List[Dict[str, Any]] = []

    def _consume_interleaved(
        items: List[Tuple[float, int, Dict[str, Any]]], start_with_light: bool
    ) -> None:
        heavy_first = sorted(items, key=lambda item: (-item[0], item[1]))
        light_first = sorted(items, key=lambda item: (item[0], item[1]))
        candidate_lists = (
            [light_first, heavy_first]
            if start_with_light
            else [heavy_first, light_first]
        )
        while True:
            progressed = False
            for candidate_list in candidate_lists:
                selected: Optional[Tuple[float, int, Dict[str, Any]]] = None
                for item in candidate_list:
                    key = item[2]["cache_key"]
                    if key in used_keys:
                        continue
                    selected = item
                    break
                if selected is None:
                    continue
                used_keys.add(selected[2]["cache_key"])
                reordered.append(selected[2])
                progressed = True
            if not progressed:
                break

    _consume_interleaved(parallel_safe, start_with_light=True)
    _consume_interleaved(parallel_unsafe, start_with_light=not parallel_safe)
    return reordered


def _adaptive_parallel_job_eligible(
    job: Dict[str, Any], dynamic_cfg: Dict[str, Any]
) -> bool:
    allowed_scenarios = {
        str(value).upper()
        for value in dynamic_cfg.get(
            "adaptive_parallel_scenarios", ["FOLDING_ONLY", "PROTEIN_ONLY"]
        )
    }
    scenario = str(job.get("scenario") or "PROTEIN_LIGAND").upper()
    if allowed_scenarios and scenario not in allowed_scenarios:
        return False
    scenario_limits = dynamic_cfg.get("adaptive_parallel_scenario_limits") or {}
    scenario_limit_cfg = scenario_limits.get(scenario) or {}
    max_estimated_vram_gb = float(
        scenario_limit_cfg.get(
            "max_job_vram_gb",
            dynamic_cfg.get("adaptive_parallel_max_job_vram_gb", 18.0),
        )
    )
    if float(job.get("estimated_vram_gb") or 0.0) > max_estimated_vram_gb:
        return False
    max_total_length = int(
        scenario_limit_cfg.get(
            "max_total_length",
            dynamic_cfg.get("adaptive_parallel_max_total_length", 700),
        )
    )
    return _boltz_job_total_sequence_length(job) <= max_total_length


def _resolve_boltz_service_predict_limit(dynamic_cfg: Dict[str, Any]) -> int:
    return max(
        1,
        int(
            dynamic_cfg.get(
                "service_max_concurrent_predicts", dynamic_cfg.get("max_workers", 1)
            )
        ),
    )


def _estimate_active_predict_vram_gb(
    active_jobs: List[Dict[str, Any]], dynamic_cfg: Dict[str, Any]
) -> float:
    if not active_jobs:
        return 0.0
    predict_limit = _resolve_boltz_service_predict_limit(dynamic_cfg)
    estimates = sorted(
        (float(job.get("estimated_vram_gb") or 0.0) for job in active_jobs),
        reverse=True,
    )
    return sum(estimates[:predict_limit])


def _resolve_adaptive_worker_target(
    boltz_client: Boltz2Client,
    dynamic_cfg: Dict[str, Any],
    logger: logging.Logger,
    current_target: int,
    active_jobs: List[Dict[str, Any]],
    next_job: Optional[Dict[str, Any]],
    telemetry_history: deque,
    low_util_streak: int,
    last_scale_at: float,
    partial_idle_since: Optional[float],
) -> Tuple[int, int, float, Optional[float]]:
    if not dynamic_cfg.get("adaptive_enabled", True):
        return (current_target, 0, last_scale_at, None)
    max_workers = int(dynamic_cfg.get("max_workers", 1))
    if max_workers <= 1:
        return (1, 0, last_scale_at, None)
    snapshot = _query_gpu_runtime_snapshot()
    if snapshot is not None:
        telemetry_history.append(snapshot)
    if not telemetry_history:
        return (1, 0, last_scale_at, partial_idle_since)
    if not active_jobs:
        return (1, 0, last_scale_at, None)
    avg_utilization = sum(
        (item["utilization_gpu"] for item in telemetry_history)
    ) / len(telemetry_history)
    peak_memory_fraction = max(
        (item["memory_utilization"] for item in telemetry_history)
    )
    recent_window = list(telemetry_history)[-min(2, len(telemetry_history)) :]
    recent_peak_utilization = max((item["utilization_gpu"] for item in recent_window))
    low_util_threshold = float(dynamic_cfg.get("scale_up_utilization_below", 35.0))
    low_memory_threshold = float(
        dynamic_cfg.get("scale_up_memory_utilization_below", 0.55)
    )
    high_memory_threshold = float(
        dynamic_cfg.get("scale_down_memory_utilization_above", 0.8)
    )
    gpu_phase_scale_up_threshold = float(
        dynamic_cfg.get("scale_up_after_gpu_activity_above", 25.0)
    )
    service_predict_limit = _resolve_boltz_service_predict_limit(dynamic_cfg)
    default_scale_up_target = 6 if service_predict_limit == 1 else 2
    scale_up_target = max(
        2,
        min(
            max_workers,
            int(dynamic_cfg.get("scale_up_target_workers", default_scale_up_target)),
        ),
    )
    required_streak = int(dynamic_cfg.get("scale_up_low_util_streak", 3))
    cooldown_seconds = float(dynamic_cfg.get("scale_cooldown_seconds", 45.0))
    min_parallel_hold_seconds = float(
        dynamic_cfg.get("min_parallel_hold_seconds", cooldown_seconds)
    )
    tail_refill_grace_seconds = float(
        dynamic_cfg.get("tail_refill_grace_seconds", 30.0)
    )
    now = time.monotonic()
    active_jobs_parallel_safe = all(
        (_adaptive_parallel_job_eligible(job, dynamic_cfg) for job in active_jobs)
    )
    next_job_parallel_safe = next_job is not None and _adaptive_parallel_job_eligible(
        next_job, dynamic_cfg
    )
    healthy = boltz_client.health_check()
    active_parallel_slots_full = len(active_jobs) >= current_target
    partial_idle = current_target > 1 and len(active_jobs) < current_target
    if partial_idle:
        partial_idle_since = partial_idle_since or now
    else:
        partial_idle_since = None
    partial_idle_elapsed = (
        0.0 if partial_idle_since is None else max(0.0, now - partial_idle_since)
    )
    insufficient_parallel_safe_work = False
    if partial_idle:
        lacks_safe_refill = not active_jobs_parallel_safe or not next_job_parallel_safe
        if lacks_safe_refill and partial_idle_elapsed >= tail_refill_grace_seconds:
            insufficient_parallel_safe_work = True
    should_scale_down = current_target > 1 and (
        not healthy
        or peak_memory_fraction >= high_memory_threshold
        or insufficient_parallel_safe_work
    )
    if should_scale_down:
        logger.info(
            "Boltz adaptive concurrency: scale down to 1 (avg_gpu_util=%.1f%% recent_peak_gpu_util=%.1f%% peak_mem=%.0f%% healthy=%s active=%d/%d partial_idle_for=%.1fs)",
            avg_utilization,
            recent_peak_utilization,
            peak_memory_fraction * 100.0,
            healthy,
            len(active_jobs),
            current_target,
            partial_idle_elapsed,
        )
        return (1, 0, now, None)
    enough_samples = len(telemetry_history) >= max(2, required_streak)
    scale_up_allowed = (
        current_target < max_workers
        and healthy
        and enough_samples
        and active_jobs_parallel_safe
        and next_job_parallel_safe
        and active_parallel_slots_full
        and (now - last_scale_at >= cooldown_seconds)
    )
    if (
        scale_up_allowed
        and avg_utilization <= low_util_threshold
        and (peak_memory_fraction <= low_memory_threshold)
    ):
        low_util_streak += 1
    else:
        low_util_streak = 0
    if (
        scale_up_allowed
        and peak_memory_fraction <= low_memory_threshold
        and (recent_peak_utilization >= gpu_phase_scale_up_threshold)
    ):
        new_target = scale_up_target
        logger.info(
            "Boltz adaptive concurrency: scale up to %d after GPU activity (recent_peak_gpu_util=%.1f%% peak_mem=%.0f%%)",
            new_target,
            recent_peak_utilization,
            peak_memory_fraction * 100.0,
        )
        return (new_target, 0, now, None)
    return (current_target, low_util_streak, last_scale_at, partial_idle_since)


def _resolve_dynamic_boltz_worker_limit(
    boltz_cfg: Dict[str, Any], pending_jobs: int
) -> int:
    dynamic_cfg = boltz_cfg.get("dynamic_concurrency") or {}
    if not dynamic_cfg.get("enabled", True):
        return max(1, int(boltz_cfg.get("workers", 1)))
    max_workers = int(
        dynamic_cfg.get("max_workers", max(1, int(boltz_cfg.get("workers", 1))))
    )
    return max(1, min(max_workers, pending_jobs))


def _build_boltz_request_metadata(
    job: Dict[str, Any],
    *,
    evo_model: str,
    boltz_artifacts_root: Path,
    boltz_artifact_record_root: Path,
    dispatch_submitted_at: Optional[str] = None,
    request_purpose: Optional[str] = None,
) -> Dict[str, Any]:
    metadata: Dict[str, Any] = {
        "mutation_id": job["mutation_id"],
        "gene_name": job["gene_name"],
        "evo_model": evo_model,
        "scenario": job["scenario"],
        "seq_hash": job["seq_hash"],
        "artifacts_dir": str(boltz_artifacts_root),
        "artifact_record_root": str(boltz_artifact_record_root),
        "fanout_mutation_ids": list(job.get("mutation_ids") or [job["mutation_id"]]),
        "estimated_duration_score": float(job.get("estimated_duration_score") or 0.0),
    }
    if dispatch_submitted_at:
        metadata["dispatch_submitted_at"] = dispatch_submitted_at
    if request_purpose:
        metadata["request_purpose"] = request_purpose
    return metadata


def _select_prefetch_ready_dispatch_index(
    queued_jobs: List[Dict[str, Any]], scan_window: int
) -> int:
    if not queued_jobs:
        return 0
    head = queued_jobs[0]
    if (
        not head.get("prefetch_submitted")
        or bool(head.get("prefetch_done"))
        or bool(head.get("prefetch_failed"))
    ):
        return 0
    for index in range(1, min(len(queued_jobs), max(1, scan_window))):
        candidate = queued_jobs[index]
        if bool(candidate.get("prefetch_done")) and (
            not bool(candidate.get("prefetch_failed"))
        ):
            return index
    return 0


def run_boltz_stage(
    config: Dict[str, Any],
    db: DatabaseManager,
    logger: logging.Logger,
    evo_model: str = "evo2_7b",
    limit: Optional[int] = None,
    gene_names: Optional[List[str]] = None,
    mutation_ids: Optional[Set[int]] = None,
) -> int:
    logger.info("=" * 80)
    logger.info("STAGE 2: BOLTZ-2 STRUCTURE PREDICTION  [model=%s]", evo_model)
    logger.info("=" * 80)
    logger.info(
        "Boltz evidence policy: core=%s | support=%s | auxiliary=%s",
        ", ".join(DatabaseManager.CORE_JUDGMENT_FIELDS[1:]),
        ", ".join(DatabaseManager.SPECIALIZED_SUPPORT_FIELDS),
        ", ".join(DatabaseManager.AUXILIARY_ARCHIVE_FIELDS[:-1]),
    )
    boltz_cfg = config["pipeline"]["boltz"]
    gene_settings = build_runtime_gene_settings(config)
    boltz_client = Boltz2Client(
        base_url=config["models"]["boltz2"]["base_url"],
        timeout_seconds=config["models"]["boltz2"]["timeout_seconds"],
        max_retries=config["models"]["boltz2"]["max_retries"],
    )
    configured_batch_size = int(boltz_cfg.get("batch_size", 32))
    processed = 0
    failed = 0
    dynamic_cfg = boltz_cfg.get("dynamic_concurrency") or {}
    service_cfg = boltz_cfg.get("service") or {}
    service_prepare_limit = max(
        1, int(service_cfg.get("max_concurrent_prepares", 1) or 1)
    )
    checkpoint_interval_seconds = int(
        config.get("pipeline", {}).get("checkpoint_sync_interval_seconds", 300)
    )
    last_checkpoint_sync = time.monotonic()
    partner_cache: Dict[Tuple[str, ...], List[str]] = {}
    export_path = Path(config["pipeline"]["export_csv_path"])
    boltz_artifacts_root = Path(
        os.getenv(
            "BOLTZ_ARTIFACTS_DIR", str(_default_boltz_artifacts_root(export_path))
        )
    )
    boltz_artifact_record_root = Path(
        os.getenv("BOLTZ_ARTIFACT_RECORD_ROOT", str(boltz_artifacts_root))
    )
    if boltz_artifacts_root != export_path.parent / "boltz_artifacts":
        logger.info(
            "Boltz artifacts default redirected to FAST storage: export_parent=%s artifacts_root=%s",
            export_path.parent,
            boltz_artifacts_root,
        )
    cache_hits = 0
    deduplicated_mutations = 0
    batch_count = 0
    prefetch_submitted_total = 0
    prefetch_completed_total = 0
    prefetch_failed_total = 0
    target_filter_sql, target_filter_params = build_target_filter_sql(
        gene_names, mutation_ids
    )
    if gene_names or mutation_ids:
        logger.info(
            "Target filter active for Boltz: genes=%s mutation_ids=%s",
            ",".join(gene_names or ["*"]),
            len(mutation_ids or set()),
        )
    while True:
        remaining_budget = None if limit is None else max(0, limit - processed - failed)
        if remaining_budget == 0:
            logger.info(
                "Boltz limit reached: attempted=%d (processed=%d, failed=%d)",
                processed + failed,
                processed,
                failed,
            )
            break
        batch_limit = (
            configured_batch_size
            if remaining_budget is None
            else min(configured_batch_size, remaining_budget)
        )
        ready_rows = db.conn.execute(
            f"\n             SELECT m.id, m.gene_id, m.pos, m.ref, m.alt, m.region_type, m.region_name, m.source_sequence,\n                 m.aa_change,\n                 m.evo_delta, g.logic_type,\n                 g.name, g.sequence, g.scenario, g.dna_sequence, g.ligand_sdf_path, g.pdb_path\n            FROM mutations m\n            JOIN genes g ON m.gene_id = g.id\n             WHERE m.status = 'BOLTZ_READY' AND m.evo_model = ?\n            {target_filter_sql}\n            ORDER BY g.id, m.pos\n            ",
            (evo_model, *target_filter_params),
        ).fetchall()
        stale_lof_rows = []
        stale_regulatory_rows = []
        ranked_candidates: List[Tuple[float, Any, str, str]] = []
        skipped_candidates = 0
        for row in ready_rows:
            logic_type = str(row["logic_type"] or "POSITIVE").upper()
            if logic_type in {
                "NEGATIVE",
                "STRUCTURAL",
            } and has_loss_of_function_signature(dict(row)):
                stale_mutation = dict(row)
                stale_mutation["status"] = "LOSS_OF_FUNCTION"
                stale_mutation["final_interpretation"] = (
                    format_loss_of_function_interpretation(stale_mutation)
                )
                stale_decision = resolve_analysis_decision(stale_mutation)
                stale_lof_rows.append(
                    (
                        stale_mutation["final_interpretation"],
                        stale_decision.resistance_phenotype,
                        stale_decision.functional_state,
                        stale_decision.evidence_code,
                        int(row["id"]),
                    )
                )
                skipped_candidates += 1
                continue
            _, scenario, region_type = _resolve_runtime_boltz_context(
                gene_settings, dict(row)
            )
            if region_type == "REGULATORY" and scenario != "PROTEIN_DNA":
                stale_mutation = dict(row)
                stale_mutation["status"] = "COMPLETED"
                stale_mutation["final_interpretation"] = (
                    f"Evo-retained regulatory candidate (delta={float(row['evo_delta'] or 0.0):.2f}); no protein-DNA Boltz scenario available, kept as Evo-only result"
                )
                stale_decision = resolve_analysis_decision(stale_mutation)
                stale_regulatory_rows.append(
                    (
                        stale_mutation["final_interpretation"],
                        stale_decision.resistance_phenotype,
                        stale_decision.functional_state,
                        stale_decision.evidence_code,
                        int(row["id"]),
                    )
                )
                skipped_candidates += 1
                continue
            score, reason, bucket = _score_boltz_candidate(
                config, gene_settings, dict(row), float(row["evo_delta"] or 0.0)
            )
            if score is None:
                skipped_candidates += 1
                continue
            ranked_candidates.append((score, row, reason, bucket))
        queue_sweep_changed = False
        if stale_lof_rows:
            db.conn.executemany(
                "\n                UPDATE mutations\n                SET status = 'LOSS_OF_FUNCTION', final_interpretation = ?,\n                    resistance_phenotype = ?, functional_state = ?, evidence_code = ?,\n                    last_error = NULL, updated_at = CURRENT_TIMESTAMP\n                WHERE id = ?\n                ",
                stale_lof_rows,
            )
            queue_sweep_changed = True
        if stale_regulatory_rows:
            db.conn.executemany(
                "\n                UPDATE mutations\n                SET status = 'COMPLETED', final_interpretation = ?,\n                    resistance_phenotype = ?, functional_state = ?, evidence_code = ?,\n                    last_error = NULL, updated_at = CURRENT_TIMESTAMP\n                WHERE id = ?\n                ",
                stale_regulatory_rows,
            )
            queue_sweep_changed = True
        if queue_sweep_changed:
            db.conn.commit()
        if stale_lof_rows:
            logger.info(
                "Boltz queue sweep: restored %d stale annotated loss-of-function BOLTZ_READY rows to LOSS_OF_FUNCTION",
                len(stale_lof_rows),
            )
        if stale_regulatory_rows:
            logger.info(
                "Boltz queue sweep: demoted %d stale regulatory BOLTZ_READY rows without a protein-DNA scenario to Evo-only COMPLETED",
                len(stale_regulatory_rows),
            )
        ranked_candidates.sort(
            key=lambda item: (
                -item[0],
                str(item[1]["name"]),
                int(item[1]["pos"]),
                int(item[1]["id"]),
            )
        )
        force_candidates = [
            item for item in ranked_candidates if item[3] == "force_boltz"
        ]
        resistance_candidates = [
            item for item in ranked_candidates if item[3] == "resistance"
        ]
        exploration_candidates = [
            item for item in ranked_candidates if item[3] == "exploration"
        ]
        target_total = batch_limit
        force_quota = min(len(force_candidates), max(0, int(target_total * 0.2)))
        resistance_quota = min(
            len(resistance_candidates), max(1, int(target_total * 0.6))
        )
        exploration_quota = min(
            len(exploration_candidates),
            max(0, target_total - force_quota - resistance_quota),
        )
        selected_items = _interleave_boltz_candidates(
            force_candidates[:force_quota],
            resistance_candidates[:resistance_quota],
            exploration_candidates[:exploration_quota],
            target_total,
        )
        selected_ids = {int(item[1]["id"]) for item in selected_items}
        if len(selected_items) < target_total:
            for item in ranked_candidates:
                item_id = int(item[1]["id"])
                if item_id in selected_ids:
                    continue
                selected_items.append(item)
                selected_ids.add(item_id)
                if len(selected_items) >= target_total:
                    break
        selected = [row for _, row, _, _ in selected_items]
        if not selected:
            if processed == 0 and failed == 0:
                logger.info("No mutations ready for Boltz. Stage 2 complete.")
            else:
                logger.info(
                    "No additional mutations ready for Boltz after %d batch(es).",
                    batch_count,
                )
            break
        batch_count += 1
        logger.info(
            "Boltz batch %d: found %d BOLTZ_READY candidates; selected top %d by dynamic priority (%d ineligible for Boltz)",
            batch_count,
            len(ready_rows),
            len(selected),
            skipped_candidates,
        )
        preview = ", ".join(
            (
                f"{row['name']}#{row['id']} score={score:.1f}"
                for score, row, _, _ in selected_items[: min(5, len(selected_items))]
            )
        )
        logger.info("Top Boltz candidates: %s", preview)
        logger.info(
            "Boltz priority mix: force_boltz=%d resistance=%d exploration=%d",
            sum((1 for item in selected_items if item[3] == "force_boltz")),
            sum((1 for item in selected_items if item[3] == "resistance")),
            sum((1 for item in selected_items if item[3] == "exploration")),
        )
        logger.info(
            "Boltz selection order pattern: 3 resistance -> 1 force_boltz -> 1 exploration (with fallback fill)"
        )
        prepared_jobs: List[Dict[str, Any]] = []
        for row in selected:
            runtime = gene_settings.get(str(row["name"]), {})
            scenario = str(
                runtime.get("scenario") or row["scenario"] or "PROTEIN_LIGAND"
            ).upper()
            region_type = str(row["region_type"] or "CDS").upper()
            try:
                if region_type == "REGULATORY":
                    if scenario != "PROTEIN_DNA":
                        raise ValueError(
                            f"Regulatory mutation for {row['name']} requires PROTEIN_DNA scenario, got {scenario}"
                        )
                    protein_sequence = translate_coding_sequence(
                        str(row["sequence"] or "")
                    )
                    dna_sequence = apply_single_base_mutation(
                        str(row["source_sequence"] or ""),
                        int(row["pos"]),
                        row["ref"],
                        str(row["alt"] or ""),
                    )
                else:
                    mutant_coding_sequence = apply_single_base_mutation(
                        str(row["source_sequence"] or row["sequence"] or ""),
                        int(row["pos"]),
                        row["ref"],
                        str(row["alt"] or ""),
                    )
                    protein_sequence = translate_coding_sequence(mutant_coding_sequence)
                    dna_sequence = (
                        str(runtime.get("dna_sequence") or row["dna_sequence"] or "")
                        or None
                    )
                partner_names = [
                    str(item) for item in runtime.get("complex_partners") or [] if item
                ]
                partner_sequences = (
                    load_partner_sequences(db, partner_names, partner_cache)
                    if partner_names
                    else []
                )
                ligand_sdf = row["ligand_sdf_path"]
                seq_hash, scenario_key, ligand_hash = _boltz_cache_key(
                    protein_sequence, scenario, ligand_sdf
                )
            except Exception as exc:
                logger.error(f"Mutation {row['id']} Boltz prep failed: {exc}")
                db.conn.execute(
                    "UPDATE mutations SET status = 'FAILED', resistance_phenotype = 'UNKNOWN', functional_state = 'UNKNOWN', evidence_code = 'ANALYSIS_FAILED', last_error = ?, retry_count = retry_count + 1, updated_at = CURRENT_TIMESTAMP WHERE id = ?",
                    (str(exc), row["id"]),
                )
                failed += 1
                last_checkpoint_sync = maybe_checkpoint_sync(
                    db,
                    logger,
                    "boltz",
                    processed,
                    failed,
                    last_checkpoint_sync,
                    checkpoint_interval_seconds,
                )
                continue
            prepared_job = {
                "mutation_id": int(row["id"]),
                "gene_name": str(row["name"]),
                "scenario": scenario,
                "scenario_key": scenario_key,
                "protein_sequence": protein_sequence,
                "protein_pdb": row["pdb_path"],
                "ligand_sdf": ligand_sdf,
                "dna_sequence": dna_sequence,
                "partner_sequences": partner_sequences,
                "seq_hash": seq_hash,
                "ligand_hash": ligand_hash,
                "cache_key": (seq_hash, scenario_key, ligand_hash),
            }
            prepared_job["estimated_vram_gb"] = _estimate_boltz_job_vram_gb(
                prepared_job, dynamic_cfg
            )
            prepared_jobs.append(prepared_job)
        if not prepared_jobs:
            continue
        prefetched_cache = _boltz_cache_prefetch(
            db, [job["cache_key"] for job in prepared_jobs]
        )
        unique_jobs: Dict[Tuple[str, str, str], Dict[str, Any]] = {}
        for job in prepared_jobs:
            cached = prefetched_cache.get(job["cache_key"])
            if cached is not None:
                cache_hits += 1
                _write_boltz_metrics_to_mutation(db, job["mutation_id"], cached)
                processed += 1
                last_checkpoint_sync = maybe_checkpoint_sync(
                    db,
                    logger,
                    "boltz",
                    processed,
                    failed,
                    last_checkpoint_sync,
                    checkpoint_interval_seconds,
                )
                if processed <= 3:
                    logger.info(
                        "Boltz cache applied to mutation %s: core=%s",
                        job["mutation_id"],
                        _format_metric_summary(
                            cached,
                            ("ptm", "iptm", "plddt", "stability_score", "clash_score"),
                        ),
                    )
                continue
            group = unique_jobs.get(job["cache_key"])
            if group is None:
                job["mutation_ids"] = [job["mutation_id"]]
                unique_jobs[job["cache_key"]] = job
                continue
            group["mutation_ids"].append(job["mutation_id"])
            deduplicated_mutations += 1
        pending_jobs = _rebalance_pending_boltz_jobs(
            list(unique_jobs.values()), dynamic_cfg
        )
        worker_limit = _resolve_dynamic_boltz_worker_limit(boltz_cfg, len(pending_jobs))
        dispatch_burst_limit = (
            worker_limit
            if not dynamic_cfg.get("adaptive_enabled", True)
            else max(
                1,
                int(
                    dynamic_cfg.get("dispatch_burst_limit", service_prepare_limit)
                    or service_prepare_limit
                ),
            )
        )
        initial_dispatch_target = (
            worker_limit
            if not dynamic_cfg.get("adaptive_enabled", True)
            else max(
                1,
                int(
                    dynamic_cfg.get(
                        "initial_dispatch_limit", min(3, service_prepare_limit)
                    )
                    or min(3, service_prepare_limit)
                ),
            )
        )
        initial_dispatch_limit = (
            worker_limit
            if not dynamic_cfg.get("adaptive_enabled", True)
            else max(
                1,
                min(
                    worker_limit,
                    dispatch_burst_limit,
                    service_prepare_limit,
                    initial_dispatch_target,
                ),
            )
        )
        msa_prefetch_cfg = service_cfg.get("msa_prefetch") or {}
        msa_prefetch_enabled = bool(msa_prefetch_cfg.get("enabled", True))
        default_prefetch_window = max(
            initial_dispatch_limit, max(4, service_prepare_limit * 2)
        )
        msa_prefetch_window = (
            0
            if not msa_prefetch_enabled
            else max(
                1,
                int(
                    msa_prefetch_cfg.get("window", default_prefetch_window)
                    or default_prefetch_window
                ),
            )
        )
        default_prefetch_workers = min(service_prepare_limit, 2)
        msa_prefetch_workers = (
            0
            if not msa_prefetch_enabled
            else max(
                1,
                min(
                    service_prepare_limit,
                    int(
                        msa_prefetch_cfg.get("max_workers", default_prefetch_workers)
                        or default_prefetch_workers
                    ),
                ),
            )
        )
        msa_prefetch_dispatch_ready_window = (
            0
            if not msa_prefetch_enabled
            else max(
                1,
                int(
                    msa_prefetch_cfg.get("dispatch_ready_window", msa_prefetch_window)
                    or msa_prefetch_window
                ),
            )
        )
        logger.info(
            "Boltz batch %d dispatch: prepared=%d cache_hits=%d unique_gpu_jobs=%d deduped=%d worker_cap=%d initial_dispatch_limit=%d burst_limit=%d",
            batch_count,
            len(prepared_jobs),
            sum((1 for job in prepared_jobs if job["cache_key"] in prefetched_cache)),
            len(pending_jobs),
            sum((max(0, len(job.get("mutation_ids", [])) - 1) for job in pending_jobs)),
            worker_limit,
            initial_dispatch_limit,
            dispatch_burst_limit,
        )
        logger.info(
            "Boltz dispatch order preview: %s",
            ", ".join(
                (
                    f"{job['gene_name']}#{job['mutation_id']} dur={float(job.get('estimated_duration_score') or 0.0):.1f} vram={float(job['estimated_vram_gb']):.1f}"
                    for job in pending_jobs[: min(6, len(pending_jobs))]
                )
            ),
        )
        logger.info(
            "Boltz MSA prefetch: enabled=%s window=%d workers=%d ready_window=%d",
            msa_prefetch_enabled,
            msa_prefetch_window,
            msa_prefetch_workers,
            msa_prefetch_dispatch_ready_window,
        )
        if not pending_jobs:
            continue
        batch_prefetch_submitted = 0
        batch_prefetch_completed = 0
        batch_prefetch_failed = 0
        prefetch_executor: Optional[ThreadPoolExecutor] = None
        try:
            if msa_prefetch_enabled and pending_jobs:
                prefetch_executor = ThreadPoolExecutor(
                    max_workers=min(msa_prefetch_workers, len(pending_jobs))
                )
            with ThreadPoolExecutor(max_workers=worker_limit) as executor:
                queued_jobs: List[Dict[str, Any]] = list(pending_jobs)
                active_futures: Dict[Any, Dict[str, Any]] = {}
                active_prefetch_futures: Dict[Any, Dict[str, Any]] = {}
                telemetry_history: deque = deque(
                    maxlen=max(1, int(dynamic_cfg.get("adaptive_telemetry_window", 6)))
                )
                adaptive_target_workers = 1
                adaptive_low_util_streak = 0
                adaptive_last_scale_at = 0.0
                adaptive_partial_idle_since: Optional[float] = None
                adaptive_poll_interval_seconds = float(
                    dynamic_cfg.get("adaptive_poll_interval_seconds", 5.0)
                )
                service_predict_limit = _resolve_boltz_service_predict_limit(
                    dynamic_cfg
                )

                def _drain_prefetch_futures(*, wait_all: bool = False) -> None:
                    nonlocal batch_prefetch_completed, batch_prefetch_failed
                    futures = list(active_prefetch_futures.keys())
                    if wait_all:
                        done_prefetch = futures
                    else:
                        done_prefetch = [future for future in futures if future.done()]
                    for future in done_prefetch:
                        job = active_prefetch_futures.pop(future)
                        try:
                            result = future.result()
                            job["prefetch_done"] = True
                            job["prefetch_failed"] = False
                            job["prefetch_result"] = result
                            batch_prefetch_completed += 1
                            if batch_prefetch_completed <= 3:
                                logger.info(
                                    "Boltz MSA prefetch ready: mutation=%s gene=%s reused_processed=%s pending=%d",
                                    job["mutation_id"],
                                    job["gene_name"],
                                    bool(result.get("reused_processed")),
                                    int(result.get("pending_records") or 0),
                                )
                        except Exception as exc:
                            job["prefetch_done"] = True
                            job["prefetch_failed"] = True
                            job["prefetch_error"] = str(exc)
                            batch_prefetch_failed += 1
                            logger.warning(
                                "Boltz MSA prefetch failed: mutation=%s gene=%s error=%s",
                                job["mutation_id"],
                                job["gene_name"],
                                exc,
                            )

                def _maybe_submit_prefetches() -> None:
                    nonlocal batch_prefetch_submitted
                    _drain_prefetch_futures(wait_all=False)
                    if prefetch_executor is None:
                        return
                    ahead_limit = max(initial_dispatch_limit, msa_prefetch_window)
                    candidate_jobs = [
                        job
                        for job in queued_jobs[: min(len(queued_jobs), ahead_limit)]
                        if not job.get("prefetch_submitted")
                        and (not job.get("prefetch_done"))
                    ]
                    while (
                        candidate_jobs
                        and len(active_prefetch_futures) < msa_prefetch_workers
                    ):
                        job = candidate_jobs.pop(0)
                        future = prefetch_executor.submit(
                            boltz_client.prefetch,
                            scenario=job["scenario"],
                            protein_sequence=job["protein_sequence"],
                            protein_pdb=job["protein_pdb"],
                            ligand_sdf=job["ligand_sdf"],
                            dna_sequence=job["dna_sequence"],
                            partner_sequences=job["partner_sequences"],
                            metadata=_build_boltz_request_metadata(
                                job,
                                evo_model=evo_model,
                                boltz_artifacts_root=boltz_artifacts_root,
                                boltz_artifact_record_root=boltz_artifact_record_root,
                                request_purpose="msa_prefetch",
                            ),
                        )
                        job["prefetch_submitted"] = True
                        active_prefetch_futures[future] = job
                        batch_prefetch_submitted += 1

                _maybe_submit_prefetches()
                while queued_jobs or active_futures:
                    _drain_prefetch_futures(wait_all=False)
                    _maybe_submit_prefetches()
                    current_budget_gb = _current_boltz_memory_budget_gb(dynamic_cfg)
                    next_job = queued_jobs[0] if queued_jobs else None
                    (
                        adaptive_target_workers,
                        adaptive_low_util_streak,
                        adaptive_last_scale_at,
                        adaptive_partial_idle_since,
                    ) = _resolve_adaptive_worker_target(
                        boltz_client,
                        dynamic_cfg,
                        logger,
                        min(worker_limit, adaptive_target_workers),
                        list(active_futures.values()),
                        next_job,
                        telemetry_history,
                        adaptive_low_util_streak,
                        adaptive_last_scale_at,
                        adaptive_partial_idle_since,
                    )
                    dispatch_limit = min(worker_limit, adaptive_target_workers)
                    if dynamic_cfg.get("adaptive_enabled", True):
                        observed_gpu_activity = any(
                            (item["utilization_gpu"] > 0 for item in telemetry_history)
                        )
                        if not observed_gpu_activity:
                            dispatch_limit = max(dispatch_limit, initial_dispatch_limit)
                    burst_submitted = 0
                    while (
                        queued_jobs
                        and len(active_futures) < dispatch_limit
                        and (burst_submitted < dispatch_burst_limit)
                    ):
                        dispatch_index = (
                            _select_prefetch_ready_dispatch_index(
                                queued_jobs, msa_prefetch_dispatch_ready_window
                            )
                            if msa_prefetch_enabled
                            else 0
                        )
                        job = queued_jobs[dispatch_index]
                        estimated_vram_gb = float(job["estimated_vram_gb"])
                        active_jobs_snapshot = list(active_futures.values())
                        active_predict_vram_gb = _estimate_active_predict_vram_gb(
                            active_jobs_snapshot, dynamic_cfg
                        )
                        if (
                            len(active_jobs_snapshot) < service_predict_limit
                            and active_predict_vram_gb + estimated_vram_gb
                            > current_budget_gb
                        ):
                            break
                        dispatch_started_at = datetime.now(timezone.utc).isoformat()
                        future = executor.submit(
                            boltz_client.predict,
                            scenario=job["scenario"],
                            protein_sequence=job["protein_sequence"],
                            protein_pdb=job["protein_pdb"],
                            ligand_sdf=job["ligand_sdf"],
                            dna_sequence=job["dna_sequence"],
                            partner_sequences=job["partner_sequences"],
                            metadata=_build_boltz_request_metadata(
                                job,
                                evo_model=evo_model,
                                boltz_artifacts_root=boltz_artifacts_root,
                                boltz_artifact_record_root=boltz_artifact_record_root,
                                dispatch_submitted_at=dispatch_started_at,
                            ),
                        )
                        job["dispatch_submitted_at"] = dispatch_started_at
                        active_futures[future] = queued_jobs.pop(dispatch_index)
                        burst_submitted += 1
                        _drain_prefetch_futures(wait_all=False)
                        _maybe_submit_prefetches()
                    if not active_futures and queued_jobs:
                        dispatch_index = (
                            _select_prefetch_ready_dispatch_index(
                                queued_jobs, msa_prefetch_dispatch_ready_window
                            )
                            if msa_prefetch_enabled
                            else 0
                        )
                        job = queued_jobs[dispatch_index]
                        dispatch_started_at = datetime.now(timezone.utc).isoformat()
                        future = executor.submit(
                            boltz_client.predict,
                            scenario=job["scenario"],
                            protein_sequence=job["protein_sequence"],
                            protein_pdb=job["protein_pdb"],
                            ligand_sdf=job["ligand_sdf"],
                            dna_sequence=job["dna_sequence"],
                            partner_sequences=job["partner_sequences"],
                            metadata=_build_boltz_request_metadata(
                                job,
                                evo_model=evo_model,
                                boltz_artifacts_root=boltz_artifacts_root,
                                boltz_artifact_record_root=boltz_artifact_record_root,
                                dispatch_submitted_at=dispatch_started_at,
                            ),
                        )
                        job["dispatch_submitted_at"] = dispatch_started_at
                        active_futures[future] = queued_jobs.pop(dispatch_index)
                        _drain_prefetch_futures(wait_all=False)
                        _maybe_submit_prefetches()
                    done, _ = wait(
                        list(active_futures.keys()),
                        timeout=max(0.5, adaptive_poll_interval_seconds),
                        return_when=FIRST_COMPLETED,
                    )
                    if not done:
                        _drain_prefetch_futures(wait_all=False)
                        _maybe_submit_prefetches()
                        continue
                    for future in done:
                        job = active_futures.pop(future)
                        mutation_ids = list(
                            job.get("mutation_ids") or [job["mutation_id"]]
                        )
                        try:
                            metrics = future.result()
                            dispatch_submitted_at = job.get("dispatch_submitted_at")
                            if dispatch_submitted_at:
                                try:
                                    elapsed = (
                                        datetime.now(timezone.utc)
                                        - datetime.fromisoformat(
                                            str(dispatch_submitted_at)
                                        )
                                    ).total_seconds()
                                    logger.info(
                                        "Boltz job end-to-end latency: mutation=%s gene=%s elapsed=%.1fs dur_score=%.1f",
                                        mutation_ids[0],
                                        job["gene_name"],
                                        elapsed,
                                        float(
                                            job.get("estimated_duration_score") or 0.0
                                        ),
                                    )
                                except ValueError:
                                    pass
                            if _boltz_metrics_empty(metrics):
                                raise ValueError(
                                    "Boltz returned no usable metrics; refusing to mark mutation as COMPLETED"
                                )
                            _boltz_cache_store(
                                db,
                                job["seq_hash"],
                                job["scenario_key"],
                                job["ligand_hash"],
                                metrics,
                            )
                            for mutation_id in mutation_ids:
                                _write_boltz_metrics_to_mutation(
                                    db, mutation_id, metrics
                                )
                                processed += 1
                            last_checkpoint_sync = maybe_checkpoint_sync(
                                db,
                                logger,
                                "boltz",
                                processed,
                                failed,
                                last_checkpoint_sync,
                                checkpoint_interval_seconds,
                            )
                            logger.debug(
                                "Boltz raw auxiliary evidence for %s: %s",
                                mutation_ids,
                                _format_metric_summary(
                                    metrics,
                                    (
                                        "binding_affinity",
                                        "affinity_pred_value",
                                        "pair_energy",
                                        "complex_energy",
                                    ),
                                ),
                            )
                            if processed <= 3:
                                logger.info(
                                    "Mutation %s Boltz completed: core=%s",
                                    mutation_ids[0],
                                    _format_metric_summary(
                                        metrics,
                                        (
                                            "ptm",
                                            "iptm",
                                            "plddt",
                                            "stability_score",
                                            "clash_score",
                                        ),
                                    ),
                                )
                            if len(mutation_ids) > 1:
                                logger.info(
                                    "Boltz dedupe fan-out applied: representative=%s copies=%d seq_hash=%s",
                                    mutation_ids[0],
                                    len(mutation_ids),
                                    job["seq_hash"][:8],
                                )
                            if processed % 10 == 0:
                                logger.info(f"Processed {processed} structures...")
                        except Exception as exc:
                            logger.error(
                                f"Mutations {mutation_ids} Boltz simulation failed: {exc}"
                            )
                            for mutation_id in mutation_ids:
                                db.conn.execute(
                                    "UPDATE mutations SET status = 'FAILED', resistance_phenotype = 'UNKNOWN', functional_state = 'UNKNOWN', evidence_code = 'ANALYSIS_FAILED', last_error = ?, retry_count = retry_count + 1, updated_at = CURRENT_TIMESTAMP WHERE id = ?",
                                    (str(exc), mutation_id),
                                )
                                failed += 1
                            last_checkpoint_sync = maybe_checkpoint_sync(
                                db,
                                logger,
                                "boltz",
                                processed,
                                failed,
                                last_checkpoint_sync,
                                checkpoint_interval_seconds,
                            )
                _drain_prefetch_futures(wait_all=True)
        finally:
            if prefetch_executor is not None:
                prefetch_executor.shutdown(wait=True)
        prefetch_submitted_total += batch_prefetch_submitted
        prefetch_completed_total += batch_prefetch_completed
        prefetch_failed_total += batch_prefetch_failed
        if msa_prefetch_enabled:
            logger.info(
                "Boltz batch %d MSA prefetch summary: submitted=%d completed=%d failed=%d",
                batch_count,
                batch_prefetch_submitted,
                batch_prefetch_completed,
                batch_prefetch_failed,
            )
    last_checkpoint_sync = maybe_checkpoint_sync(
        db,
        logger,
        "boltz",
        processed,
        failed,
        last_checkpoint_sync,
        checkpoint_interval_seconds,
        force=True,
    )
    logger.info(
        "Stage 2 complete: %d processed (%d cache hits, %d deduped), %d failed across %d batch(es) [msa_prefetch submitted=%d completed=%d failed=%d]",
        processed,
        cache_hits,
        deduplicated_mutations,
        failed,
        batch_count,
        prefetch_submitted_total,
        prefetch_completed_total,
        prefetch_failed_total,
    )
    append_batch_metrics(
        Path(config["pipeline"]["metrics_path"]),
        {
            "stage": "boltz",
            "processed": processed,
            "failed": failed,
            "msa_prefetch_submitted": prefetch_submitted_total,
            "msa_prefetch_completed": prefetch_completed_total,
            "msa_prefetch_failed": prefetch_failed_total,
        },
    )
    return 0 if failed == 0 else 1


def run_analyze_stage(
    config: Dict[str, Any], db: DatabaseManager, logger: logging.Logger
) -> int:
    logger.info("=" * 80)
    logger.info("STAGE 3: ANALYSIS & EXPORT")
    logger.info("=" * 80)
    db.clear_synergy_events()
    pathway_rules = config.get("pathway_rules", [])
    synergy_rules = config.get("synergy_rules", [])
    analyzer = PathwayAnalyzer(db, pathway_rules, synergy_rules)
    completed = db.conn.execute(
        "\n         SELECT m.id, m.gene_id, g.name AS gene_name, g.logic_type, m.region_type, m.aa_change, m.ref, m.alt,\n                 m.evo_delta, m.boltz_iptm, m.boltz_plddt, m.boltz_affinity_pred_value,\n             m.boltz_ptm, m.boltz_stability_score, m.synergy_score, m.synergy_tags,\n                 m.final_interpretation, m.status,\n                 g.sequence AS gene_sequence, g.scenario, g.ligand_sdf_path,\n             m.resistance_phenotype, m.functional_state, m.evidence_code\n        FROM mutations m\n        JOIN genes g ON m.gene_id = g.id\n           WHERE m.status IN ('COMPLETED', 'LOSS_OF_FUNCTION', 'LETHAL_SKIP')\n        "
    ).fetchall()
    logger.info(f"Found {len(completed)} analyzable mutations for analysis")
    batch = [dict(row) for row in completed]
    batch = analyzer.apply_epistasis_corrections(batch)
    batch = analyzer.apply_multi_gene_synergy(batch)
    for mut in batch:
        db.conn.execute(
            "\n            UPDATE mutations\n            SET status = ?, final_interpretation = ?, synergy_score = ?, synergy_tags = ?,\n                resistance_phenotype = ?, functional_state = ?, evidence_code = ?,\n                updated_at = CURRENT_TIMESTAMP\n            WHERE id = ?\n            ",
            (
                mut.get("status"),
                mut.get("final_interpretation"),
                mut.get("synergy_score", 0.0),
                mut.get("synergy_tags"),
                mut.get("resistance_phenotype"),
                mut.get("functional_state"),
                mut.get("evidence_code"),
                mut["id"],
            ),
        )
    db.conn.commit()
    export_path = Path(config["pipeline"]["export_csv_path"])
    export_path.parent.mkdir(parents=True, exist_ok=True)
    results = db.conn.execute(
        "\n        SELECT g.name, g.logic_type, g.scenario, m.aa_change, m.region_type, m.region_name, m.pos, m.ref, m.alt,\n               m.evo_model, m.evo_delta,\n             m.boltz_binding_affinity, m.boltz_affinity_pred_value, m.boltz_pair_energy,\n             m.boltz_stability_score, m.boltz_ptm,\n             m.boltz_iptm, m.boltz_plddt, m.boltz_complex_energy,\n             m.synergy_score, m.synergy_tags,\n               m.status, m.final_interpretation,\n               m.resistance_phenotype, m.functional_state, m.evidence_code\n        FROM mutations m\n        JOIN genes g ON m.gene_id = g.id\n        WHERE m.status IN ('COMPLETED', 'LOSS_OF_FUNCTION', 'RE_SENSITIZED')\n        ORDER BY g.name, m.evo_model, m.id\n        "
    ).fetchall()
    _write_analysis_csv(export_path, results, logger, "Results")
    analysis_snapshot_path_value = str(
        config["pipeline"].get("analysis_snapshot_csv_path", "")
    ).strip()
    if analysis_snapshot_path_value:
        analysis_snapshot_path = Path(analysis_snapshot_path_value)
        if analysis_snapshot_path != export_path:
            _write_analysis_csv(
                analysis_snapshot_path,
                results,
                logger,
                "Git-friendly analysis snapshot",
            )
    db.sync_to_share()
    logger.info("Stage 3 complete")
    return 0


def main() -> int:
    parser = ArgumentParser(
        description="AMR-Hunter: Antibiotic Resistance Prediction Pipeline"
    )
    parser.add_argument(
        "--stage",
        choices=["evo", "boltz", "analyze", "all"],
        default="all",
        help="Pipeline stage to execute (default: all)",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path(__file__).parent / "config" / "config.yaml",
        help="Path to config.yaml",
    )
    parser.add_argument(
        "--evo-model",
        default=None,
        dest="evo_model",
        help="Evo2 model name to use for scoring (e.g. evo2_7b, evo2_20b). Defaults to config models.evo2.model_name or env EVO2_MODEL_NAME, then 'evo2_7b'.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Limit the number of mutations processed per stage (useful for test runs).",
    )
    parser.add_argument(
        "--gene",
        action="append",
        dest="genes",
        default=None,
        help="Restrict Evo/Boltz stages to a target gene. Repeat for multiple genes.",
    )
    parser.add_argument(
        "--mutation-id-file",
        type=Path,
        default=None,
        help="Restrict Evo/Boltz stages to mutation ids listed in a text file.",
    )
    args = parser.parse_args()
    config = load_config(args.config)
    import os

    evo_model: str = (
        args.evo_model
        or os.getenv("EVO2_MODEL_NAME")
        or str(
            (config.get("models", {}) or {}).get("evo2", {}).get("model_name")
            or "evo2_7b"
        )
    )
    limit: Optional[int] = args.limit
    gene_names = normalize_gene_filter(args.genes)
    try:
        mutation_ids = load_mutation_id_filter(args.mutation_id_file)
    except ValueError as exc:
        print(f"Invalid --mutation-id-file: {exc}", file=sys.stderr)
        return 2
    logger = setup_logger(Path(config["pipeline"]["log_path"]))
    logger.info(
        "Active evo model: %s%s", evo_model, f"  [limit={limit}]" if limit else ""
    )
    if gene_names or mutation_ids:
        logger.info(
            "Targeted run filter: genes=%s mutation_ids=%s",
            ",".join(gene_names or ["*"]),
            len(mutation_ids or set()),
        )
    db = DatabaseManager(
        local_db_path=config["paths"]["local_db_path"],
        share_db_path=config["paths"]["share_db_path"],
    )
    try:
        if args.stage in ("evo", "all"):
            if (
                run_evo_stage(
                    config,
                    db,
                    logger,
                    evo_model=evo_model,
                    limit=limit,
                    gene_names=gene_names,
                    mutation_ids=mutation_ids,
                )
                != 0
            ):
                logger.error("Evo stage failed")
                return 1
        if args.stage in ("boltz", "all"):
            if (
                run_boltz_stage(
                    config,
                    db,
                    logger,
                    evo_model=evo_model,
                    limit=limit,
                    gene_names=gene_names,
                    mutation_ids=mutation_ids,
                )
                != 0
            ):
                logger.error("Boltz stage failed")
                return 1
        if args.stage in ("analyze", "all"):
            if run_analyze_stage(config, db, logger) != 0:
                logger.error("Analyze stage failed")
                return 1
        logger.info("Pipeline execution successful")
        return 0
    except Exception as exc:
        logger.exception(f"Fatal error: {exc}")
        return 1
    finally:
        db.close()


if __name__ == "__main__":
    sys.exit(main())
