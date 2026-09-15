from __future__ import annotations
import re
from dataclasses import dataclass
from typing import Any, Mapping, Optional, Set

BDQ_LOF_NEUTRAL_GENES: Set[str] = {"glpK", "Rv1979c"}
BDQ_RV0678_NEUTRAL_FS_POSITIONS: Set[int] = {64, 71, 78, 10, 44, 124, 140, 164}


@dataclass(frozen=True)
class ConsequenceDecision:
    resistance_phenotype: str
    functional_state: str
    evidence_code: str


def normalize_consequence_label(value: Any) -> str:
    return str(value or "").strip().lower().replace("_", "-").replace(" ", "-")


def is_loss_of_function_consequence(value: Any) -> bool:
    text = str(value or "").strip()
    if not text:
        return False
    normalized = normalize_consequence_label(text)
    if normalized in {
        "lof",
        "loss-of-function",
        "loss-of-function-variant",
        "frameshift",
        "frameshift-variant",
        "start-lost",
        "start-lost-variant",
        "stop-gained",
        "stop-gain",
        "nonsense",
    }:
        return True
    lower = text.lower()
    if "frameshift" in lower or "fs" in lower:
        return True
    if normalized.startswith("p.met1") or normalized.endswith("?"):
        return True
    if re.fullmatch("[Mm]1[A-Za-z]", text):
        return True
    return text.endswith("*") or text.endswith("Ter")


def is_bdq_lof_resistant_gene(gene: str, logic_type: str) -> bool:
    gene_name = str(gene or "").strip()
    normalized_logic = str(logic_type or "").strip().upper()
    if not gene_name or not normalized_logic:
        return False
    if gene_name in BDQ_LOF_NEUTRAL_GENES:
        return False
    return normalized_logic == "NEGATIVE"


def infer_loss_of_function_decision(
    logic_type: Any, consequence: Any = "loss_of_function"
) -> ConsequenceDecision:
    normalized_logic = str(logic_type or "POSITIVE").strip().upper()
    if normalized_logic == "NEGATIVE":
        return ConsequenceDecision(
            "RESISTANT", "LOSS_OF_FUNCTION", "NEGATIVE_LOF_SIGNATURE"
        )
    if normalized_logic == "STRUCTURAL":
        return ConsequenceDecision(
            "SENSITIVE", "LOSS_OF_FUNCTION", "STRUCTURAL_LOF_SIGNATURE"
        )
    if normalized_logic == "POSITIVE":
        return ConsequenceDecision(
            "SENSITIVE", "LOSS_OF_FUNCTION", "POSITIVE_TARGET_COLLAPSE"
        )
    return ConsequenceDecision("UNKNOWN", "LOSS_OF_FUNCTION", "UNKNOWN_LOGIC_LOF")


def _parse_fs_position(mutation: Mapping[str, Any]) -> Optional[int]:
    for field in ("mutation", "aa_change"):
        value = str(mutation.get(field) or "").strip()
        match = re.fullmatch("p\\.?[A-Za-z*?]+(\\d+)fs", value, re.IGNORECASE)
        if match:
            return int(match.group(1))
    return None


def resolve_consequence_level_decision(
    mutation: Mapping[str, Any], drug: str = "", ignore_rv0678_neutral_fs: bool = False
) -> Optional[ConsequenceDecision]:
    candidates = (
        mutation.get("consequence"),
        mutation.get("mutation_class"),
        mutation.get("aa_change"),
        mutation.get("mutation"),
    )
    if not any((is_loss_of_function_consequence(value) for value in candidates)):
        return None
    gene = str(mutation.get("gene") or "")
    logic_type = str(mutation.get("logic_type") or "")
    if drug == "Bedaquiline" and gene and (logic_type == "NEGATIVE"):
        if not is_bdq_lof_resistant_gene(gene, logic_type):
            return None
        if gene == "Rv0678" and (not ignore_rv0678_neutral_fs):
            fs_pos = _parse_fs_position(mutation)
            if fs_pos is not None and fs_pos in BDQ_RV0678_NEUTRAL_FS_POSITIONS:
                return None
    return infer_loss_of_function_decision(logic_type, "loss_of_function")
