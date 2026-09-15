from __future__ import annotations
import hashlib
from dataclasses import dataclass
import logging
from pathlib import Path
import re
from typing import Any, Dict, List, Mapping, Optional, Set, Tuple
from Bio.Seq import Seq
from core.consequence import (
    BDQ_LOF_NEUTRAL_GENES,
    infer_loss_of_function_decision,
    is_loss_of_function_consequence,
    normalize_consequence_label,
    resolve_consequence_level_decision,
)
from core.database import DatabaseManager

TARGET_GENES = {
    "atpE",
    "Rv0678",
    "mmpS5",
    "mmpL5",
    "pepQ",
    "lpqB",
    "mtrA",
    "mtrB",
    "Rv1979c",
    "glpK",
}


@dataclass(frozen=True)
class AnalysisDecision:
    resistance_phenotype: str
    functional_state: str
    evidence_code: str


def _coerce_float(value: Any) -> Optional[float]:
    try:
        return None if value is None else float(value)
    except (TypeError, ValueError):
        return None


def _translate_coding_sequence(sequence: str) -> Optional[str]:
    normalized = (sequence or "").upper().replace("U", "T")
    if not normalized or len(normalized) % 3 != 0:
        return None
    translated = str(Seq(normalized).translate(table="11", to_stop=False))
    if "*" in translated:
        translated = translated.split("*", 1)[0]
    return translated or None


def _compute_affinity_shift(mutation: Mapping[str, Any]) -> Optional[float]:
    affinity_pred_value = _coerce_float(mutation.get("boltz_affinity_pred_value"))
    wt_affinity_pred_value = _coerce_float(mutation.get("wt_boltz_affinity_pred_value"))
    if affinity_pred_value is not None and wt_affinity_pred_value is not None:
        return affinity_pred_value - wt_affinity_pred_value
    interpretation = str(mutation.get("final_interpretation") or "")
    shift_match = re.search("affinity_shift=([-\\d.]+)", interpretation)
    if shift_match:
        try:
            return float(shift_match.group(1))
        except (ValueError, TypeError):
            pass
    mt_match = re.search("mt_affinity=([-\\d.]+)", interpretation)
    wt_match = re.search("wt_affinity=([-\\d.]+)", interpretation)
    if mt_match and wt_match:
        try:
            return float(mt_match.group(1)) - float(wt_match.group(1))
        except (ValueError, TypeError):
            pass
    return None


def is_synonymous_cds_variant(mutation: Mapping[str, Any]) -> bool:
    if str(mutation.get("region_type") or "CDS").upper() != "CDS":
        return False
    aa_change = str(mutation.get("aa_change") or "").strip()
    if not aa_change:
        return False
    simple_match = re.fullmatch("([A-Z*?])(\\d+)([A-Z*?])", aa_change)
    if simple_match is not None:
        return simple_match.group(1) == simple_match.group(3)
    hgvs_match = re.fullmatch("p\\.([A-Za-z*?]+)(\\d+)([A-Za-z*?]+)", aa_change)
    if hgvs_match is not None:
        return hgvs_match.group(1).lower() == hgvs_match.group(3).lower()
    return False


def format_loss_of_function_interpretation(mutation: Mapping[str, Any]) -> str:
    logic_type = str(mutation.get("logic_type") or "POSITIVE").upper()
    aa_change = mutation.get("aa_change") or "<unknown>"
    if logic_type == "NEGATIVE":
        return f"Resistant by annotated truncating/frameshift/start-loss signature: aa_change={aa_change}"
    if logic_type == "STRUCTURAL":
        return f"SENSITIVE: annotated truncating/frameshift/start-loss signature likely collapses protein function: aa_change={aa_change}"
    return f"SENSITIVE: annotated truncating/frameshift/start-loss signature likely destroys the drug target: aa_change={aa_change}"


def _normalized_synergy_tags(mutation: Mapping[str, Any]) -> Set[str]:
    raw_tags = str(mutation.get("synergy_tags") or "")
    return {tag for tag in raw_tags.split("|") if tag}


def _has_boltz_core(mutation: Mapping[str, Any]) -> bool:
    return any(
        (
            mutation.get(field) is not None
            for field in ("boltz_ptm", "boltz_iptm", "boltz_plddt")
        )
    )


def _is_supported_regulatory_scenario(mutation: Mapping[str, Any]) -> bool:
    scenario = str(mutation.get("scenario") or "").strip().upper()
    return scenario == "PROTEIN_DNA"


def _normalized_scenario(mutation: Mapping[str, Any]) -> str:
    return str(mutation.get("scenario") or "").strip().upper()


def _is_rv0678_regulatory_binding_escape(
    mutation: Mapping[str, Any], core_score: Optional[float], plddt: Optional[float]
) -> bool:
    gene_name = str(mutation.get("gene_name") or mutation.get("name") or "")
    region_type = str(mutation.get("region_type") or "CDS").strip().upper()
    if gene_name != "Rv0678" or region_type != "REGULATORY":
        return False
    if _normalized_scenario(mutation) != "PROTEIN_DNA":
        return False
    if core_score is None or plddt is None:
        return False
    return (
        core_score <= PathwayAnalyzer.RV0678_DNA_BINDING_MAX_CORE
        and plddt >= PathwayAnalyzer.RV0678_DNA_BINDING_MIN_PLDDT
    )


def _is_efflux_complex_interface_disrupted(
    mutation: Mapping[str, Any], core_score: Optional[float], plddt: Optional[float]
) -> bool:
    gene_name = str(mutation.get("gene_name") or mutation.get("name") or "")
    region_type = str(mutation.get("region_type") or "CDS").strip().upper()
    if (
        gene_name not in PathwayAnalyzer.EFFLUX_COMPLEX_INTERFACE_GENES
        or region_type != "CDS"
    ):
        return False
    if _normalized_scenario(mutation) != "COMPLEX":
        return False
    if core_score is None or plddt is None:
        return False
    return (
        core_score <= PathwayAnalyzer.EFFLUX_COMPLEX_INTERFACE_MAX_CORE
        and plddt >= PathwayAnalyzer.EFFLUX_COMPLEX_INTERFACE_MIN_PLDDT
    )


def resolve_analysis_decision(mutation: Mapping[str, Any]) -> AnalysisDecision:
    _gene_role = str(mutation.get("gene_role") or "").strip().lower()
    _head_gene_name = str(
        mutation.get("gene_name") or mutation.get("name") or ""
    ).strip()
    if _gene_role == "exploratory" or _head_gene_name in BDQ_LOF_NEUTRAL_GENES:
        _head_status = str(mutation.get("status") or "").strip().upper()
        if _head_status == "LOSS_OF_FUNCTION":
            _exploratory_functional = "LOSS_OF_FUNCTION"
        elif _head_status == "LETHAL_SKIP":
            _exploratory_functional = "LETHAL"
        elif _head_status == "COMPLETED":
            _exploratory_functional = "INTACT"
        else:
            _exploratory_functional = "UNKNOWN"
        return AnalysisDecision("UNKNOWN", _exploratory_functional, "EXPLORATORY_ONLY")
    stored_resistance = str(mutation.get("resistance_phenotype") or "").strip().upper()
    stored_functional = str(mutation.get("functional_state") or "").strip().upper()
    stored_evidence = str(mutation.get("evidence_code") or "").strip().upper()
    _stale_evidence_codes = {
        "ANALYSIS_FAILED",
        "PENDING_ANALYSIS",
        "QUEUED_FOR_BOLTZ",
        "ATPE_AFFINITY_REQUIRED",
        "UNRESOLVED",
    }
    _current_status = str(mutation.get("status") or "").strip().upper()
    _current_logic_type = str(mutation.get("logic_type") or "POSITIVE").strip().upper()
    if stored_resistance and stored_functional and stored_evidence:
        _evidence_logic_stale = _current_logic_type != "NEGATIVE" and (
            stored_evidence.startswith("NEGATIVE_")
            or stored_evidence == "MTRAB_DEREGULATION"
        )
        if not _evidence_logic_stale and (
            stored_evidence not in _stale_evidence_codes
            or _current_status not in ("COMPLETED", "LOSS_OF_FUNCTION", "RE_SENSITIZED")
        ):
            return AnalysisDecision(
                stored_resistance, stored_functional, stored_evidence
            )
    status = str(mutation.get("status") or "").strip().upper()
    logic_type = str(mutation.get("logic_type") or "POSITIVE").strip().upper()
    region_type = str(mutation.get("region_type") or "CDS").strip().upper()
    gene_name = str(mutation.get("gene_name") or mutation.get("name") or "").strip()
    interpretation = str(mutation.get("final_interpretation") or "").strip().lower()
    synergy_tags = _normalized_synergy_tags(mutation)
    evo_delta = _coerce_float(mutation.get("evo_delta"))
    ptm = _coerce_float(mutation.get("boltz_ptm"))
    iptm = _coerce_float(mutation.get("boltz_iptm"))
    plddt = _coerce_float(mutation.get("boltz_plddt"))
    affinity_shift = _coerce_float(mutation.get("boltz_affinity_shift"))
    if affinity_shift is None:
        affinity_shift = _compute_affinity_shift(mutation)
    core_values = [value for value in (ptm, iptm) if value is not None]
    core_score = min(core_values) if core_values else None
    if is_synonymous_cds_variant(mutation):
        return AnalysisDecision("SENSITIVE", "INTACT", "SYNONYMOUS_NO_PROTEIN_CHANGE")
    if status == "FAILED":
        return AnalysisDecision("UNKNOWN", "UNKNOWN", "ANALYSIS_FAILED")
    if status in {"PENDING", "EVO_DONE"}:
        return AnalysisDecision("UNKNOWN", "UNKNOWN", "PENDING_ANALYSIS")
    if status == "BOLTZ_READY":
        return AnalysisDecision("UNKNOWN", "UNKNOWN", "QUEUED_FOR_BOLTZ")
    if status == "RE_SENSITIZED":
        return AnalysisDecision("SENSITIVE", "RE_SENSITIZED", "PATHWAY_RE_SENSITIZED")
    if "efflux_pump_collapse" in synergy_tags or "PUMP_COLLAPSE" in synergy_tags:
        functional_state = (
            "LOSS_OF_FUNCTION" if status == "LOSS_OF_FUNCTION" else "INTACT"
        )
        return AnalysisDecision("SENSITIVE", functional_state, "EFFLUX_PUMP_COLLAPSE")
    if (
        "atp_synthase_compensation" in synergy_tags
        or "ATP_COMPENSATION" in synergy_tags
    ) and (not has_loss_of_function_signature(dict(mutation))):
        return AnalysisDecision("RESISTANT", "INTACT", "ATP_SYNTHASE_COMPENSATION")
    if status == "LETHAL_SKIP":
        return AnalysisDecision("SENSITIVE", "LETHAL", "EVO_LETHAL")
    if status == "LOSS_OF_FUNCTION":
        if logic_type == "NEGATIVE":
            evidence_code = (
                "NEGATIVE_LOF_SIGNATURE"
                if has_loss_of_function_signature(dict(mutation))
                else "NEGATIVE_LOF_REVIEWED"
            )
            return AnalysisDecision("RESISTANT", "LOSS_OF_FUNCTION", evidence_code)
        if logic_type == "STRUCTURAL":
            evidence_code = (
                "STRUCTURAL_LOF_SIGNATURE"
                if has_loss_of_function_signature(dict(mutation))
                else "STRUCTURAL_COLLAPSE"
            )
            return AnalysisDecision("SENSITIVE", "LOSS_OF_FUNCTION", evidence_code)
        return AnalysisDecision(
            "SENSITIVE", "LOSS_OF_FUNCTION", "POSITIVE_TARGET_COLLAPSE"
        )
    if status == "COMPLETED":
        if region_type == "REGULATORY":
            if not _is_supported_regulatory_scenario(mutation):
                return AnalysisDecision(
                    "UNKNOWN", "UNKNOWN", "REGULATORY_SCENARIO_UNSUPPORTED"
                )
            if not _has_boltz_core(mutation):
                return AnalysisDecision("UNKNOWN", "UNKNOWN", "REGULATORY_EVO_ONLY")
        if _is_rv0678_regulatory_binding_escape(mutation, core_score, plddt):
            return AnalysisDecision(
                "RESISTANT", "LOSS_OF_FUNCTION", "RV0678_DNA_BINDING_ESCAPE"
            )
        if _is_efflux_complex_interface_disrupted(mutation, core_score, plddt):
            return AnalysisDecision(
                "SENSITIVE", "LOSS_OF_FUNCTION", "EFFLUX_COMPLEX_INTERFACE_DISRUPTED"
            )
        if logic_type == "NEGATIVE":
            return AnalysisDecision("SENSITIVE", "INTACT", "NEGATIVE_FUNCTION_RETAINED")
        if logic_type == "STRUCTURAL":
            return AnalysisDecision(
                "SENSITIVE", "INTACT", "STRUCTURAL_FUNCTION_RETAINED"
            )
        if gene_name not in TARGET_GENES:
            return AnalysisDecision("UNKNOWN", "INTACT", "NON_TARGET_GENE")
        if (
            gene_name == "atpE"
            and (not has_loss_of_function_signature(dict(mutation)))
            and (evo_delta is not None)
            and (evo_delta <= PathwayAnalyzer.ATPE_RESISTANCE_EVO_THRESHOLD)
            and (core_score is not None)
            and (plddt is not None)
            and (core_score >= PathwayAnalyzer.ATPE_RESISTANCE_MIN_CORE)
            and (plddt >= PathwayAnalyzer.ATPE_RESISTANCE_MIN_PLDDT)
        ):
            if affinity_shift is not None:
                if affinity_shift >= PathwayAnalyzer.ATPE_RESISTANCE_MIN_AFFINITY_SHIFT:
                    return AnalysisDecision(
                        "RESISTANT", "INTACT", "ATPE_TARGET_ESCAPE_AFFINITY_SHIFT"
                    )
                return AnalysisDecision(
                    "UNKNOWN", "INTACT", "ATPE_AFFINITY_NO_ESCAPE_SHIFT"
                )
            return AnalysisDecision("UNKNOWN", "INTACT", "ATPE_AFFINITY_REQUIRED")
        return AnalysisDecision(
            "UNKNOWN", "INTACT", "POSITIVE_REVIEWED_NO_ESCAPE_EVIDENCE"
        )
    if interpretation.startswith("resistant"):
        return AnalysisDecision("RESISTANT", "UNKNOWN", "INTERPRETATION_PREFIX")
    if interpretation.startswith("sensitive"):
        return AnalysisDecision("SENSITIVE", "UNKNOWN", "INTERPRETATION_PREFIX")
    return AnalysisDecision("UNKNOWN", "UNKNOWN", "UNRESOLVED")


def attach_analysis_decision(mutation: Dict[str, Any]) -> AnalysisDecision:
    decision = resolve_analysis_decision(mutation)
    mutation["resistance_phenotype"] = decision.resistance_phenotype
    mutation["functional_state"] = decision.functional_state
    mutation["evidence_code"] = decision.evidence_code
    return decision


def derive_resistance_label(
    status: Any, logic_type: Any, final_interpretation: Any
) -> str:
    return resolve_analysis_decision(
        {
            "status": status,
            "logic_type": logic_type,
            "final_interpretation": final_interpretation,
        }
    ).resistance_phenotype


def has_loss_of_function_signature(mutation: Dict[str, Any]) -> bool:
    if str(mutation.get("region_type") or "CDS").upper() != "CDS":
        return False
    if any(
        (
            is_loss_of_function_consequence(mutation.get(field))
            for field in ("aa_change", "consequence", "mutation_class", "mutation")
        )
    ):
        return True
    ref = str(mutation.get("ref") or "")
    alt = str(mutation.get("alt") or "")
    return bool(ref and alt and (len(ref) != len(alt)))


class PathwayAnalyzer:
    EFFLUX_STABILITY_SCALE = 10.0
    ACTIVE_ANALYSIS_STATUSES = {"COMPLETED"}
    RV0678_COLLAPSE_MAX_CORE = 0.75
    RV0678_COLLAPSE_MAX_PLDDT = 0.83
    RV0678_DNA_BINDING_MAX_CORE = 0.78
    RV0678_DNA_BINDING_MIN_PLDDT = 0.84
    ATPE_RESISTANCE_MIN_CORE = 0.92
    ATPE_RESISTANCE_MIN_PLDDT = 0.94
    ATPE_RESISTANCE_EVO_THRESHOLD = -8.0
    ATPE_RESISTANCE_MIN_AFFINITY_SHIFT = 0.08
    EFFLUX_COMPLEX_INTERFACE_GENES = frozenset({"mmpS5", "mmpL5"})
    EFFLUX_COMPLEX_INTERFACE_MAX_CORE = 0.7
    EFFLUX_COMPLEX_INTERFACE_MIN_PLDDT = 0.85

    def __init__(
        self,
        db: DatabaseManager,
        pathway_rules: Optional[List[Dict[str, Any]]] = None,
        synergy_rules: Optional[List[Dict[str, Any]]] = None,
    ) -> None:
        self.db = db
        self.pathway_rules = pathway_rules or []
        self.synergy_rules = synergy_rules or []
        self.logger = logging.getLogger("amr_hunter.analyzer")
        self._pathway_index = self._build_pathway_index()
        self._wt_boltz_metrics_cache: Dict[
            Tuple[str, str, str], Optional[Dict[str, Optional[float]]]
        ] = {}
        self._ligand_hash_cache: Dict[str, str] = {}

    def _build_pathway_index(self) -> Dict[str, Dict[str, Any]]:
        index: Dict[str, Dict[str, Any]] = {}
        for rule in self.pathway_rules:
            pathway_group = rule.get("pathway_group")
            if not pathway_group:
                continue
            index[pathway_group] = {
                "regulator": rule.get("regulator"),
                "effectors": set(rule.get("effectors", [])),
                "epistasis_logic": rule.get("epistasis_logic"),
            }
        return index

    def build_rv0678_efflux_network(self) -> Dict[str, Any]:
        network = {
            "regulator": "Rv0678",
            "effectors": {"mmpS5", "mmpL5"},
            "resistant_statuses": {"LOSS_OF_FUNCTION"},
            "collapsed_statuses": {"LOSS_OF_FUNCTION", "LETHAL_SKIP"},
            "flip_to": "SENSITIVE",
        }
        for rule in self.pathway_rules:
            regulator = rule.get("regulator")
            effectors = set(rule.get("effectors", []))
            if regulator == "Rv0678" and {"mmpS5", "mmpL5"}.issubset(effectors):
                return {
                    **network,
                    "effectors": effectors,
                    "flip_to": rule.get("correction_state", "SENSITIVE"),
                }
        return network

    def check_epistasis(
        self, gene_name: str, status: str, logic_type: str
    ) -> Optional[Tuple[str, str]]:
        pathway_group = None
        for pg, info in self._pathway_index.items():
            if gene_name == info.get("regulator") or gene_name in info.get(
                "effectors", set()
            ):
                pathway_group = pg
                break
        if not pathway_group:
            return None
        pathway_info = self._pathway_index[pathway_group]
        regulator = pathway_info["regulator"]
        effectors = pathway_info["effectors"]
        if gene_name == regulator and status not in self.ACTIVE_ANALYSIS_STATUSES:
            return self._check_regulator_epistasis(regulator, effectors, logic_type)
        if gene_name in effectors:
            return self._check_effector_epistasis(regulator, effectors, logic_type)
        return None

    def _check_regulator_epistasis(
        self, regulator: str, effectors: Set[str], logic_type: str
    ) -> Optional[Tuple[str, str]]:
        effector_statuses = self.db.conn.execute(
            "\n                 SELECT g.name, COUNT(CASE WHEN m.status IN ('LOSS_OF_FUNCTION', 'LETHAL_SKIP') THEN 1 END) as loss_count,\n                     COUNT(CASE WHEN m.status IN ('COMPLETED') THEN 1 END) as completed_count\n            FROM genes g\n            LEFT JOIN mutations m ON m.gene_id = g.id\n            WHERE g.name IN ({})\n            GROUP BY g.id\n            ".format(
                ",".join((f"'{e}'" for e in effectors))
            )
        ).fetchall()
        if not effector_statuses:
            return None
        disrupted = sum((row["loss_count"] for row in effector_statuses))
        total = sum(
            (row["loss_count"] + row["completed_count"] for row in effector_statuses)
        )
        if total > 0 and disrupted / total > 0.5:
            return (
                "RE_SENSITIZED",
                f"Regulator {regulator} and effectors {effectors} co-disrupted, pathway non-functional",
            )
        return None

    def _check_effector_epistasis(
        self, regulator: str, effectors: Set[str], logic_type: str
    ) -> Optional[Tuple[str, str]]:
        regulator_data = self.db.conn.execute(
            "\n                 SELECT COUNT(CASE WHEN m.status IN ('LOSS_OF_FUNCTION', 'LETHAL_SKIP') THEN 1 END) as loss_count,\n                     COUNT(CASE WHEN m.status IN ('COMPLETED') THEN 1 END) as completed_count\n            FROM genes g\n            LEFT JOIN mutations m ON m.gene_id = g.id\n            WHERE g.name = ?\n            ",
            (regulator,),
        ).fetchone()
        if not regulator_data:
            return None
        regulator_intact = (
            regulator_data["completed_count"] > 0 and regulator_data["loss_count"] == 0
        )
        return None

    def apply_epistasis_corrections(
        self, mutations_batch: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        for mutation in mutations_batch:
            self._restore_stale_loss_of_function(mutation)
            self._attach_affinity_context(mutation)
            self._apply_gene_specific_completed_inference(mutation)
            self._ensure_completed_interpretation(mutation)
            gene_name = str(mutation.get("gene_name") or mutation.get("name") or "")
            status = str(mutation.get("status") or "")
            logic_type = str(mutation.get("logic_type") or "POSITIVE")
            correction = self.check_epistasis(gene_name, status, logic_type)
            if correction:
                corrected_status, reason = correction
                mutation["status"] = corrected_status
                mutation["final_interpretation"] = reason
                self._record_pathway_correction_event(
                    regulator_mutation=mutation,
                    reason=reason,
                    corrected_status=corrected_status,
                )
                self.logger.info(
                    f"Epistasis correction: {gene_name} -> {corrected_status} [{reason}]"
                )
        self._apply_rv0678_efflux_flip(mutations_batch)
        for mutation in mutations_batch:
            attach_analysis_decision(mutation)
        return mutations_batch

    def _record_pathway_correction_event(
        self, regulator_mutation: Dict[str, Any], reason: str, corrected_status: str
    ) -> Optional[int]:
        gene_name = str(
            regulator_mutation.get("gene_name") or regulator_mutation.get("name") or ""
        )
        pathway_group = None
        for pg, info in self._pathway_index.items():
            if gene_name == info.get("regulator") or gene_name in info.get(
                "effectors", set()
            ):
                pathway_group = pg
                break
        if not pathway_group:
            return None
        pathway_info = self._pathway_index[pathway_group]
        effectors = pathway_info.get("effectors", set())
        event_id = self._record_synergy_event(
            rule_name=f"pathway_epistasis_{gene_name}",
            rule_type="REGULATOR_EFFECTOR_EPISTASIS",
            rule_source="pathway_rules",
            description=reason,
            score_modifier=0.0,
        )
        self._link_event_mutation(
            event_id,
            regulator_mutation["id"],
            "causal",
            str(regulator_mutation.get("status") or ""),
        )
        if effectors:
            placeholders = ",".join(("?" for _ in effectors))
            effector_lof_rows = self.db.conn.execute(
                f"SELECT m.id, g.name AS gene_name, m.status\n                    FROM mutations m\n                    JOIN genes g ON m.gene_id = g.id\n                    WHERE g.name IN ({placeholders})\n                      AND m.status IN ('LOSS_OF_FUNCTION', 'LETHAL_SKIP')",
                tuple(effectors),
            ).fetchall()
            for ef_row in effector_lof_rows:
                self._link_event_mutation(
                    event_id, ef_row["id"], "causal", ef_row["status"]
                )
        pathway_genes = {gene_name} | effectors
        if len(pathway_genes) > 1:
            placeholders = ",".join(("?" for _ in pathway_genes))
            other_rows = self.db.conn.execute(
                f"SELECT m.id, m.status\n                    FROM mutations m\n                    JOIN genes g ON m.gene_id = g.id\n                    WHERE g.name IN ({placeholders})\n                      AND m.status IN ('COMPLETED', 'LOSS_OF_FUNCTION', 'LETHAL_SKIP')\n                      AND m.id != ?",
                tuple(pathway_genes) + (regulator_mutation["id"],),
            ).fetchall()
            for row in other_rows:
                self._link_event_mutation(
                    event_id, row["id"], "affected", row["status"]
                )
        return event_id

    def _restore_stale_loss_of_function(self, mutation: Dict[str, Any]) -> None:
        if mutation.get("status") != "COMPLETED":
            return
        logic_type = str(mutation.get("logic_type") or "POSITIVE").upper()
        if logic_type not in {"NEGATIVE", "STRUCTURAL"}:
            return
        if not has_loss_of_function_signature(mutation):
            return
        mutation["status"] = "LOSS_OF_FUNCTION"
        mutation["final_interpretation"] = format_loss_of_function_interpretation(
            mutation
        )

    def _compute_ligand_hash(self, ligand_sdf_path: Optional[str]) -> str:
        cache_key = str(ligand_sdf_path or "")
        if cache_key in self._ligand_hash_cache:
            return self._ligand_hash_cache[cache_key]
        ligand_payload = ""
        if cache_key:
            try:
                ligand_payload = Path(cache_key).read_text(encoding="utf-8")
            except OSError:
                ligand_payload = cache_key
        ligand_hash = (
            hashlib.md5(ligand_payload.encode("utf-8")).hexdigest()
            if ligand_payload
            else ""
        )
        self._ligand_hash_cache[cache_key] = ligand_hash
        return ligand_hash

    def _lookup_wildtype_boltz_metrics(
        self, mutation: Mapping[str, Any]
    ) -> Optional[Dict[str, Optional[float]]]:
        gene_sequence = str(mutation.get("gene_sequence") or "").strip()
        scenario = str(mutation.get("scenario") or "").strip().upper()
        if not gene_sequence or not scenario:
            return None
        wt_protein = _translate_coding_sequence(gene_sequence)
        if not wt_protein:
            return None
        ligand_hash = self._compute_ligand_hash(
            str(mutation.get("ligand_sdf_path") or "").strip() or None
        )
        cache_key = (
            hashlib.md5(wt_protein.encode("utf-8")).hexdigest(),
            scenario,
            ligand_hash,
        )
        if cache_key in self._wt_boltz_metrics_cache:
            return self._wt_boltz_metrics_cache[cache_key]
        row = self.db.conn.execute(
            "\n            SELECT affinity_pred_value, iptm, ptm, plddt\n            FROM boltz_cache\n            WHERE seq_hash = ? AND scenario = ? AND ligand_hash = ?\n            ",
            cache_key,
        ).fetchone()
        metrics = None
        if row is not None:
            metrics = {
                "affinity_pred_value": _coerce_float(row["affinity_pred_value"]),
                "iptm": _coerce_float(row["iptm"]),
                "ptm": _coerce_float(row["ptm"]),
                "plddt": _coerce_float(row["plddt"]),
            }
        self._wt_boltz_metrics_cache[cache_key] = metrics
        return metrics

    def _attach_affinity_context(self, mutation: Dict[str, Any]) -> None:
        if mutation.get("status") != "COMPLETED":
            return
        if str(mutation.get("gene_name") or mutation.get("name") or "") != "atpE":
            return
        affinity_pred_value = self._coerce_float(
            mutation.get("boltz_affinity_pred_value")
        )
        if affinity_pred_value is None:
            return
        wt_metrics = self._lookup_wildtype_boltz_metrics(mutation)
        if not wt_metrics:
            return
        wt_affinity_pred_value = _coerce_float(wt_metrics.get("affinity_pred_value"))
        if wt_affinity_pred_value is None:
            return
        mutation["wt_boltz_affinity_pred_value"] = wt_affinity_pred_value
        mutation["boltz_affinity_shift"] = affinity_pred_value - wt_affinity_pred_value

    def _apply_gene_specific_completed_inference(
        self, mutation: Dict[str, Any]
    ) -> None:
        if mutation.get("status") != "COMPLETED":
            return
        gene_name = str(mutation.get("gene_name") or mutation.get("name") or "")
        region_type = str(mutation.get("region_type") or "CDS").upper()
        scenario = _normalized_scenario(mutation)
        ptm = self._coerce_float(mutation.get("boltz_ptm"))
        iptm = self._coerce_float(mutation.get("boltz_iptm"))
        plddt = self._coerce_float(mutation.get("boltz_plddt"))
        evo_delta = self._coerce_float(mutation.get("evo_delta"))
        affinity_pred_value = self._coerce_float(
            mutation.get("boltz_affinity_pred_value")
        )
        wt_affinity_pred_value = self._coerce_float(
            mutation.get("wt_boltz_affinity_pred_value")
        )
        affinity_shift = self._coerce_float(mutation.get("boltz_affinity_shift"))
        if affinity_shift is None:
            affinity_shift = _compute_affinity_shift(mutation)
            if affinity_shift is not None:
                mutation["boltz_affinity_shift"] = affinity_shift
        core_values = [value for value in (ptm, iptm) if value is not None]
        core_score = min(core_values) if core_values else None
        if core_score is None or plddt is None:
            return
        if (
            gene_name == "Rv0678"
            and core_score <= self.RV0678_COLLAPSE_MAX_CORE
            and (plddt <= self.RV0678_COLLAPSE_MAX_PLDDT)
        ):
            mutation["status"] = "LOSS_OF_FUNCTION"
            mutation["final_interpretation"] = (
                f"Resistant: Rv0678 likely loses regulatory function after Boltz review (ptm={core_score:.3f}, plddt={plddt:.3f})"
            )
            return
        if _is_rv0678_regulatory_binding_escape(mutation, core_score, plddt):
            mutation["status"] = "LOSS_OF_FUNCTION"
            mutation["final_interpretation"] = (
                f"Resistant: Rv0678 regulatory variant is predicted to disrupt protein-DNA binding while retaining overall fold confidence (ptm={core_score:.3f}, plddt={plddt:.3f}, scenario={scenario}, region={region_type})"
            )
            return
        if _is_efflux_complex_interface_disrupted(mutation, core_score, plddt):
            mutation["status"] = "LOSS_OF_FUNCTION"
            mutation["final_interpretation"] = (
                f"SENSITIVE: efflux complex interface is predicted to be disrupted despite preserved monomer confidence (gene={gene_name}, ptm={core_score:.3f}, plddt={plddt:.3f}, scenario={scenario})"
            )
            return
        if (
            gene_name == "atpE"
            and evo_delta is not None
            and (evo_delta <= self.ATPE_RESISTANCE_EVO_THRESHOLD)
            and (core_score >= self.ATPE_RESISTANCE_MIN_CORE)
            and (plddt >= self.ATPE_RESISTANCE_MIN_PLDDT)
        ):
            existing = str(mutation.get("final_interpretation") or "").strip()
            if affinity_shift is not None:
                if affinity_shift >= self.ATPE_RESISTANCE_MIN_AFFINITY_SHIFT:
                    resistant_text = f"Resistant: atpE target-altering variant remains structurally viable and shows reduced predicted binding versus WT after Boltz review (delta={evo_delta:.2f}, ptm={core_score:.3f}, plddt={plddt:.3f}, affinity_shift={affinity_shift:.3f}, mt_affinity={affinity_pred_value:.3f}, wt_affinity={wt_affinity_pred_value:.3f})"
                    mutation["final_interpretation"] = (
                        f"{resistant_text}; {existing}".strip("; ")
                    )
                    return
                mutation["final_interpretation"] = (
                    f"REVIEWED: atpE variant remains structurally viable, but predicted binding shift versus WT is below the escape threshold (delta={evo_delta:.2f}, ptm={core_score:.3f}, plddt={plddt:.3f}, affinity_shift={affinity_shift:.3f}, mt_affinity={affinity_pred_value:.3f}, wt_affinity={wt_affinity_pred_value:.3f})"
                )
                return
            reviewed_text = f"REVIEWED: atpE target-altering variant remains structurally viable after Boltz review, but affinity context versus WT is unavailable so no resistance call is assigned (delta={evo_delta:.2f}, ptm={core_score:.3f}, plddt={plddt:.3f}"
            if affinity_pred_value is not None:
                reviewed_text += (
                    f", mt_affinity={affinity_pred_value:.3f}, wt_affinity=unavailable)"
                )
            else:
                reviewed_text += ", affinity=unavailable)"
            mutation["final_interpretation"] = f"{reviewed_text}; {existing}".strip(
                "; "
            )

    @staticmethod
    def _coerce_float(value: Any) -> Optional[float]:
        return _coerce_float(value)

    def _ensure_completed_interpretation(self, mutation: Dict[str, Any]) -> None:
        if mutation.get("status") != "COMPLETED":
            return
        if mutation.get("final_interpretation"):
            return
        logic_type = str(mutation.get("logic_type") or "POSITIVE").upper()
        region_type = str(mutation.get("region_type") or "CDS").upper()
        gene_name = str(mutation.get("gene_name") or mutation.get("name") or "")
        scenario = _normalized_scenario(mutation)
        if region_type == "REGULATORY" and (
            not _is_supported_regulatory_scenario(mutation)
        ):
            mutation["final_interpretation"] = (
                "UNSUPPORTED: regulatory candidate has no protein-DNA Boltz scenario and is retained as Evo-only evidence"
            )
            return
        if (
            gene_name == "Rv0678"
            and region_type == "REGULATORY"
            and (scenario == "PROTEIN_DNA")
        ):
            mutation["final_interpretation"] = (
                "REVIEWED: Rv0678 regulatory protein-DNA candidate retained complex confidence; no calibrated DNA-binding escape threshold was crossed"
            )
            return
        if (
            gene_name in self.EFFLUX_COMPLEX_INTERFACE_GENES
            and region_type == "CDS"
            and (scenario == "COMPLEX")
        ):
            mutation["final_interpretation"] = (
                "REVIEWED: efflux complex candidate retained interface confidence; no calibrated complex-collapse threshold was crossed"
            )
            return
        if logic_type == "NEGATIVE":
            mutation["final_interpretation"] = (
                "SENSITIVE: no loss-of-function evidence after Evo/Boltz review"
            )
            return
        if logic_type == "STRUCTURAL":
            mutation["final_interpretation"] = (
                "SENSITIVE: no structural collapse evidence after Evo/Boltz review"
            )
            return
        if gene_name not in TARGET_GENES:
            mutation["final_interpretation"] = (
                "REVIEWED: gene lies outside the target panel; evidence retained without a calibrated resistance rule"
            )
            return
        mutation["final_interpretation"] = (
            "REVIEWED: Evo/Boltz review completed; no rule-based resistance label assigned"
        )

    def apply_multi_gene_synergy(
        self, mutations_batch: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        if not self.synergy_rules or not mutations_batch:
            return mutations_batch
        gene_status_map = self._build_gene_status_map(mutations_batch)
        for rule in self.synergy_rules:
            if not self._rule_matches(rule, gene_status_map):
                continue
            rule_name = str(rule.get("name") or "UNNAMED_SYNERGY")
            rule_type = str(rule.get("rule_type") or rule.get("logic") or "CO_LOF")
            score = float(
                rule.get("score", rule.get("synergy_score_modifier", 1.0)) or 1.0
            )
            interpretation = str(
                rule.get("interpretation")
                or rule.get("description")
                or f"Multi-gene synergy matched: {rule_name}"
            )
            apply_status = rule.get("apply_status")
            eligible_current_statuses = set(
                rule.get("eligible_current_statuses", ["COMPLETED"])
            )
            configured_tags = set(rule.get("tags", []))
            affected_genes = set(rule.get("genes", []))
            regulator = rule.get("regulator")
            if regulator:
                affected_genes.add(regulator)
            effectors = list(rule.get("effectors", []))
            if not effectors and rule.get("effector"):
                effectors = [rule.get("effector")]
            affected_genes.update(effectors)
            event_id = self._record_synergy_event(
                rule_name=rule_name,
                rule_type=rule_type,
                rule_source="synergy_rules",
                description=interpretation,
                score_modifier=score,
            )
            causal_map = self._build_causal_status_map(rule, affected_genes)
            for row in mutations_batch:
                gene_name = row.get("gene_name") or row.get("name")
                if affected_genes and gene_name not in affected_genes:
                    continue
                row_status = str(row.get("status") or "")
                causal_set = causal_map.get(gene_name, set())
                role = "causal" if row_status in causal_set else "affected"
                self._link_event_mutation(event_id, row["id"], role, row_status)
                tags = set(filter(None, str(row.get("synergy_tags") or "").split("|")))
                tags.add(rule_name)
                tags.update(configured_tags)
                row["synergy_tags"] = "|".join(sorted(tags))
                row["synergy_score"] = float(row.get("synergy_score") or 0.0) + score
                existing = row.get("final_interpretation") or ""
                row["final_interpretation"] = f"{existing}; {interpretation}".strip(
                    "; "
                )
                if apply_status and row.get("status") in eligible_current_statuses:
                    row["status"] = str(apply_status)
            self.logger.info(
                "Synergy rule matched: %s (event_id=%s)", rule_name, event_id
            )
        for row in mutations_batch:
            attach_analysis_decision(row)
        return mutations_batch

    def _build_causal_status_map(
        self, rule: Dict[str, Any], affected_genes: Set[str]
    ) -> Dict[str, Set[str]]:
        logic = str(rule.get("logic") or rule.get("rule_type") or "CO_LOF")
        causal_map: Dict[str, Set[str]] = {}
        if logic == "CO_LOF":
            trigger = set(
                rule.get("collapse_statuses", ["LOSS_OF_FUNCTION", "LETHAL_SKIP"])
            )
            for gene in affected_genes:
                causal_map[gene] = trigger
        elif logic == "CO_COMPLETED":
            for gene in affected_genes:
                causal_map[gene] = self.ACTIVE_ANALYSIS_STATUSES
        elif logic == "REGULATOR_EFFECTOR_DEREPRESSION":
            regulator = str(rule.get("regulator") or "")
            regulator_trigger = set(
                rule.get("regulator_statuses", ["LOSS_OF_FUNCTION"])
            )
            effector_trigger = set(rule.get("effector_active_statuses", ["COMPLETED"]))
            effectors = list(rule.get("effectors", []))
            if not effectors and rule.get("effector"):
                effectors = [str(rule.get("effector"))]
            effector_set = {str(e) for e in effectors if e}
            for gene in affected_genes:
                if gene == regulator:
                    causal_map[gene] = regulator_trigger
                elif gene in effector_set:
                    causal_map[gene] = effector_trigger
                else:
                    causal_map[gene] = set()
        return causal_map

    def _record_synergy_event(
        self,
        rule_name: str,
        rule_type: str,
        rule_source: str,
        description: str,
        score_modifier: float,
    ) -> int:
        cursor = self.db.conn.execute(
            "INSERT INTO synergy_events (rule_name, rule_type, rule_source, description, synergy_score_modifier)\n               VALUES (?, ?, ?, ?, ?)",
            (rule_name, rule_type, rule_source, description, score_modifier),
        )
        self.db.conn.commit()
        return cursor.lastrowid

    def _link_event_mutation(
        self, event_id: int, mutation_id: int, role: str, status_at_match: str
    ) -> None:
        self.db.conn.execute(
            "INSERT OR IGNORE INTO synergy_event_mutations (event_id, mutation_id, role, status_at_match)\n               VALUES (?, ?, ?, ?)",
            (event_id, mutation_id, role, status_at_match),
        )

    def _build_gene_status_map(
        self, mutations_batch: List[Dict[str, Any]]
    ) -> Dict[str, Set[str]]:
        index: Dict[str, Set[str]] = {}
        for row in mutations_batch:
            gene_name = row.get("gene_name") or row.get("name")
            status = row.get("status")
            if not gene_name or not status:
                continue
            index.setdefault(gene_name, set()).add(str(status))
        return index

    def _rule_matches(
        self, rule: Dict[str, Any], gene_status_map: Dict[str, Set[str]]
    ) -> bool:
        logic = str(rule.get("logic") or rule.get("rule_type") or "CO_LOF")
        if logic == "CO_LOF":
            genes = rule.get("genes", [])
            min_hits = int(rule.get("min_hits", 2) or 2)
            collapse_statuses = set(
                rule.get("collapse_statuses", ["LOSS_OF_FUNCTION", "LETHAL_SKIP"])
            )
            hits = sum(
                (
                    1
                    for gene in genes
                    if collapse_statuses.intersection(gene_status_map.get(gene, set()))
                )
            )
            return hits >= min_hits
        if logic == "CO_COMPLETED":
            genes = rule.get("genes", [])
            min_hits = int(rule.get("min_hits", 2) or 2)
            hits = sum(
                (
                    1
                    for gene in genes
                    if self.ACTIVE_ANALYSIS_STATUSES.intersection(
                        gene_status_map.get(gene, set())
                    )
                )
            )
            return hits >= min_hits
        if logic == "REGULATOR_EFFECTOR_DEREPRESSION":
            regulator = rule.get("regulator")
            effectors = list(rule.get("effectors", []))
            if not effectors and rule.get("effector"):
                effectors = [rule.get("effector")]
            normalized_effectors = [str(effector) for effector in effectors if effector]
            regulator_statuses = set(
                rule.get("regulator_statuses", ["LOSS_OF_FUNCTION"])
            )
            effector_active_statuses = set(
                rule.get("effector_active_statuses", ["COMPLETED"])
            )
            regulator_match = bool(
                regulator
                and regulator_statuses.intersection(
                    gene_status_map.get(regulator, set())
                )
            )
            effector_match = any(
                (
                    effector_active_statuses.intersection(
                        gene_status_map.get(effector, set())
                    )
                    for effector in normalized_effectors
                )
            )
            return regulator_match and effector_match
        return False

    def _apply_rv0678_efflux_flip(self, mutations_batch: List[Dict[str, Any]]) -> None:
        network = self.build_rv0678_efflux_network()
        regulator = network["regulator"]
        effectors = network["effectors"]
        effector_rows = [
            row
            for row in mutations_batch
            if (row.get("gene_name") or row.get("name")) in effectors
        ]
        if not effector_rows:
            return
        loss_rates = [
            self._estimate_effector_loss_rate(row, network["collapsed_statuses"])
            for row in effector_rows
        ]
        effector_loss_rate = sum(loss_rates) / len(loss_rates) if loss_rates else 0.0
        if effector_loss_rate <= 0:
            return
        event_id = self._record_synergy_event(
            rule_name="rv0678_efflux_flip",
            rule_type="REGULATOR_EFFECTOR_FLIP",
            rule_source="pathway_rules",
            description=f"Rv0678 resistance discounted by efflux collapse (loss_rate={effector_loss_rate:.2f})",
            score_modifier=-effector_loss_rate,
        )
        for row in mutations_batch:
            if (row.get("gene_name") or row.get("name")) == regulator and row.get(
                "status"
            ) in network["resistant_statuses"]:
                weighted_score = max(0.0, 1.0 - effector_loss_rate)
                row["synergy_score"] = (
                    float(row.get("synergy_score") or 0.0) - effector_loss_rate
                )
                tags = set(filter(None, str(row.get("synergy_tags") or "").split("|")))
                tags.add("EFFLUX_WEIGHTED")
                row["synergy_tags"] = "|".join(sorted(tags))
                self._link_event_mutation(
                    event_id, row["id"], "causal", str(row.get("status") or "")
                )
                if effector_loss_rate >= 0.95:
                    row["status"] = "RE_SENSITIZED"
                    row["final_interpretation"] = (
                        f"SENSITIVE: Rv0678 resistance is canceled by near-complete efflux pump collapse (weighted loss={effector_loss_rate:.2f})"
                    )
                else:
                    existing = row.get("final_interpretation") or ""
                    addition = f"Rv0678 resistance discounted by explicit efflux stability-loss evidence: final_score={weighted_score:.2f}, loss_rate={effector_loss_rate:.2f}"
                    row["final_interpretation"] = f"{existing}; {addition}".strip("; ")
        for row in effector_rows:
            gene_status = str(row.get("status") or "")
            role = (
                "causal" if gene_status in network["collapsed_statuses"] else "affected"
            )
            self._link_event_mutation(event_id, row["id"], role, gene_status)

    def _estimate_effector_loss_rate(
        self, row: Dict[str, Any], collapsed_statuses: Set[str]
    ) -> float:
        if row.get("status") in collapsed_statuses:
            return 1.0
        stability_score = row.get("boltz_stability_score")
        if stability_score is None:
            return 0.0
        try:
            numeric = float(stability_score)
        except (TypeError, ValueError):
            return 0.0
        if numeric >= 0:
            return 0.0
        return max(0.0, min(1.0, -numeric / self.EFFLUX_STABILITY_SCALE))


def apply_differential_filtering(
    mutations_batch: List[Dict[str, Any]],
    evo_cfg: Dict[str, float],
    gene_settings: Optional[Dict[str, Dict[str, Any]]] = None,
) -> List[Dict[str, Any]]:
    for mutation in mutations_batch:
        logic_type = mutation.get("logic_type", "POSITIVE")
        evo_delta = mutation.get("evo_delta")
        gene_name = str(mutation.get("name") or mutation.get("gene_name") or "")
        force_boltz = bool(
            (gene_settings or {}).get(gene_name, {}).get("force_boltz", False)
        )
        if evo_delta is None:
            continue
        lethal_threshold = evo_cfg.get("lethal_threshold", -10.0)
        collapse_threshold = evo_cfg.get("collapse_threshold", -8.0)
        if is_synonymous_cds_variant(mutation):
            mutation["status"] = "COMPLETED"
            mutation["final_interpretation"] = (
                "SENSITIVE: synonymous CDS variant leaves protein sequence unchanged in the current protein-centric model"
            )
            attach_analysis_decision(mutation)
            continue
        if logic_type in {"NEGATIVE", "STRUCTURAL"} and has_loss_of_function_signature(
            mutation
        ):
            if force_boltz:
                mutation["status"] = "BOLTZ_READY"
                mutation["final_interpretation"] = (
                    f"Force-Boltz override for annotated loss-of-function variant: aa_change={mutation.get('aa_change') or '<unknown>'}"
                )
            else:
                mutation["status"] = "LOSS_OF_FUNCTION"
                mutation["final_interpretation"] = (
                    format_loss_of_function_interpretation(mutation)
                )
            attach_analysis_decision(mutation)
            continue
        if logic_type == "POSITIVE":
            if evo_delta < lethal_threshold:
                if force_boltz:
                    mutation["status"] = "BOLTZ_READY"
                    mutation["final_interpretation"] = (
                        f"Force-Boltz override for high-value target: delta={evo_delta:.2f} < {lethal_threshold}"
                    )
                else:
                    mutation["status"] = "LETHAL_SKIP"
                    mutation["final_interpretation"] = (
                        f"Lethal mutation (delta={evo_delta:.2f} < {lethal_threshold})"
                    )
            else:
                mutation["status"] = "BOLTZ_READY"
        elif logic_type == "NEGATIVE":
            if evo_delta < collapse_threshold:
                if force_boltz:
                    mutation["status"] = "BOLTZ_READY"
                    mutation["final_interpretation"] = (
                        f"Force-Boltz override for NEGATIVE logic: delta={evo_delta:.2f} < {collapse_threshold}; run Boltz before final LOF call"
                    )
                else:
                    mutation["status"] = "LOSS_OF_FUNCTION"
                    mutation["final_interpretation"] = (
                        f"Resistant by short-circuit (NEGATIVE logic): function collapse delta={evo_delta:.2f} < {collapse_threshold}; skip Boltz-2"
                    )
            else:
                mutation["status"] = "BOLTZ_READY"
        elif logic_type == "STRUCTURAL":
            if evo_delta < lethal_threshold:
                if force_boltz:
                    mutation["status"] = "BOLTZ_READY"
                    mutation["final_interpretation"] = (
                        f"Force-Boltz override for structural target: delta={evo_delta:.2f} < {lethal_threshold}"
                    )
                else:
                    mutation["status"] = "LETHAL_SKIP"
                    mutation["final_interpretation"] = (
                        f"Lethal structural collapse (delta={evo_delta:.2f} < {lethal_threshold})"
                    )
            elif evo_delta < collapse_threshold:
                if force_boltz:
                    mutation["status"] = "BOLTZ_READY"
                    mutation["final_interpretation"] = (
                        f"Force-Boltz override for partial structural collapse: delta={evo_delta:.2f} < {collapse_threshold}"
                    )
                else:
                    mutation["status"] = "LOSS_OF_FUNCTION"
                    mutation["final_interpretation"] = (
                        f"Structural efflux collapse (delta={evo_delta:.2f} < {collapse_threshold}); potential re-sensitization context"
                    )
            else:
                mutation["status"] = "BOLTZ_READY"
        else:
            mutation["status"] = "BOLTZ_READY"
        attach_analysis_decision(mutation)
    return mutations_batch
