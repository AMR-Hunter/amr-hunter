from __future__ import annotations
import argparse
import re
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple
import pandas as pd
import yaml

AA3 = {
    "A": "Ala",
    "R": "Arg",
    "N": "Asn",
    "D": "Asp",
    "C": "Cys",
    "Q": "Gln",
    "E": "Glu",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "L": "Leu",
    "K": "Lys",
    "M": "Met",
    "F": "Phe",
    "P": "Pro",
    "S": "Ser",
    "T": "Thr",
    "W": "Trp",
    "Y": "Tyr",
    "V": "Val",
    "*": "*",
}
FINAL_STATUSES = {"COMPLETED", "LOSS_OF_FUNCTION", "RE_SENSITIZED"}
DEFAULT_CONFIG = Path(
    "/root/gpufree-share/amr_hunter_workspace/code/config/config.yaml"
)
MODELED_GENES = [
    "Rv0678",
    "atpE",
    "pepQ",
    "mtrA",
    "mtrB",
    "mmpL5",
    "mmpS5",
    "glpK",
    "Rv1979c",
]
GENERIC_NON_EXECUTABLE_MUTATIONS = {"LoF"}
SIMPLE_CDS_MUTATION_RE = re.compile("c\\.(\\d+)([ACGT]+)>([ACGT]+)")
UPSTREAM_MUTATION_RE = re.compile("c\\.\\-(\\d+)([ACGT]+)>([ACGT]+)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare AMR-Hunter predictions against WHO catalogue entries."
    )
    parser.add_argument(
        "--db",
        type=Path,
        default=Path("/root/gpufree-data/db/amr_hunter.db"),
        help="Path to SQLite database.",
    )
    parser.add_argument(
        "--xlsx",
        type=Path,
        default=Path(
            "/root/gpufree-share/amr_hunter_workspace/amr_hunter_data/reference/who/WHO-UCN-TB-2023.5-eng.xlsx"
        ),
        help="Path to WHO workbook.",
    )
    parser.add_argument(
        "--drug",
        default="Bedaquiline",
        help="Drug name in the WHO workbook to evaluate.",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG,
        help="Path to config.yaml used to resolve regulatory support windows.",
    )
    parser.add_argument(
        "--truth",
        choices=["all", "resistant", "sensitive"],
        default="all",
        help="Restrict evaluation to resistant or sensitive WHO rows.",
    )
    parser.add_argument(
        "--extra-modeled-gene",
        action="append",
        default=[],
        help="Append an optional gene to the default modeled-gene set. Repeatable.",
    )
    return parser.parse_args()


def resolve_modeled_genes(extra_genes: Iterable[str]) -> List[str]:
    modeled = list(MODELED_GENES)
    seen = set(modeled)
    for gene in extra_genes:
        normalized = str(gene or "").strip()
        if not normalized or normalized in seen:
            continue
        modeled.append(normalized)
        seen.add(normalized)
    return modeled


def truth_label(text: object) -> Optional[str]:
    if not isinstance(text, str):
        return None
    normalized = text.strip()
    if normalized.startswith(("1)", "2)")):
        return "R"
    if normalized.startswith(("4)", "5)")):
        return "S"
    return None


def is_strict_comparable_cds_variant(mutation: object, effect: object) -> bool:
    if not isinstance(mutation, str):
        return False
    mutation = mutation.strip()
    if isinstance(effect, str) and effect == "upstream_gene_variant":
        return False
    if mutation in GENERIC_NON_EXECUTABLE_MUTATIONS:
        return False
    if mutation.startswith("p."):
        return all(
            (token not in mutation for token in ("del", "ins", "dup", "fs", "?"))
        )
    return bool(SIMPLE_CDS_MUTATION_RE.fullmatch(mutation))


def normalize_window_specs(config_path: Path) -> Dict[str, List[Dict[str, Any]]]:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    default_windows = list((config.get("noncoding") or {}).get("default_windows") or [])
    window_specs: Dict[str, List[Dict[str, Any]]] = {}
    for item in config.get("target_genes") or []:
        gene = str(item.get("gene") or "").strip()
        if not gene:
            continue
        raw_windows: List[Dict[str, Any]] = []
        if bool(item.get("apply_default_windows", False)):
            raw_windows.extend(default_windows)
        raw_windows.extend(item.get("noncoding_windows") or [])
        if bool(item.get("include_upstream", False)):
            upstream_bp = max(1, int(item.get("upstream_bp", 120) or 120))
            raw_windows.append(
                {
                    "name": f"UPSTREAM_{upstream_bp}BP",
                    "mode": "upstream",
                    "length_bp": upstream_bp,
                }
            )
        normalized: List[Dict[str, Any]] = []
        for window in raw_windows:
            mode = str(window.get("mode", "upstream") or "upstream").lower()
            name = str(window.get("name") or f"{gene}_{mode}")
            spec: Dict[str, Any] = {"name": name, "mode": mode}
            if mode == "upstream":
                spec["length_bp"] = max(
                    1,
                    int(window.get("length_bp", window.get("upstream_bp", 120)) or 120),
                )
            elif mode == "around_start":
                spec["upstream_bp"] = max(0, int(window.get("upstream_bp", 80) or 80))
                spec["downstream_bp"] = max(
                    0, int(window.get("downstream_bp", 20) or 20)
                )
            else:
                continue
            normalized.append(spec)
        window_specs[gene] = normalized
    return window_specs


def regulatory_candidate_keys(
    gene: str, mutation: object, gene_window_specs: Dict[str, List[Dict[str, Any]]]
) -> Set[Tuple[str, str, int, str, str]]:
    if not isinstance(mutation, str):
        return set()
    parsed = UPSTREAM_MUTATION_RE.fullmatch(mutation.strip())
    if not parsed:
        return set()
    upstream_distance = int(parsed.group(1))
    ref = parsed.group(2)
    alt = parsed.group(3)
    candidate_keys: Set[Tuple[str, str, int, str, str]] = set()
    for spec in gene_window_specs.get(gene, []):
        mode = str(spec.get("mode") or "")
        name = str(spec.get("name") or "REGULATORY")
        if mode == "upstream":
            length_bp = int(spec.get("length_bp") or 0)
            if 1 <= upstream_distance <= length_bp:
                candidate_keys.add(
                    (gene, name, length_bp - upstream_distance + 1, ref, alt)
                )
        elif mode == "around_start":
            upstream_bp = int(spec.get("upstream_bp") or 0)
            if 1 <= upstream_distance <= upstream_bp:
                candidate_keys.add(
                    (gene, name, upstream_bp - upstream_distance + 1, ref, alt)
                )
    return candidate_keys


def classify_current_support(
    gene: object,
    mutation: object,
    effect: object,
    gene_window_specs: Dict[str, List[Dict[str, Any]]],
) -> Tuple[bool, Optional[str]]:
    if not isinstance(gene, str) or not isinstance(mutation, str):
        return (False, "missing_variant_identity")
    mutation = mutation.strip()
    if mutation in GENERIC_NON_EXECUTABLE_MUTATIONS:
        return (False, "generic_lof_category")
    if isinstance(effect, str) and effect == "upstream_gene_variant":
        if regulatory_candidate_keys(gene, mutation, gene_window_specs):
            return (True, None)
        return (False, "regulatory_out_of_config_window")
    if mutation.startswith("p."):
        return (True, None)
    if SIMPLE_CDS_MUTATION_RE.fullmatch(mutation):
        return (True, None)
    return (False, "unsupported_hgvs_form")


def aa_change_to_mutations(change: object) -> Set[str]:
    if not isinstance(change, str) or len(change) < 3:
        return set()
    if change.startswith("p."):
        if change.endswith("Ter"):
            return {change, f"{change[:-3]}*"}
        if change.endswith("*"):
            return {change, f"{change[:-1]}Ter"}
        return {change}
    wt = change[0]
    mt = change[-1]
    pos = change[1:-1]
    if wt not in AA3 or mt not in AA3 or (not pos.isdigit()):
        return set()
    mutation = f"p.{AA3[wt]}{int(pos)}{AA3[mt]}"
    if mt == "*":
        return {mutation, f"p.{AA3[wt]}{int(pos)}Ter"}
    return {mutation}


def load_who_catalog(
    xlsx_path: Path,
    drug: str,
    truth_mode: str,
    gene_window_specs: Dict[str, List[Dict[str, Any]]],
) -> pd.DataFrame:
    who = pd.read_excel(
        xlsx_path,
        sheet_name="Catalogue_master_file",
        header=2,
        usecols=["drug", "gene", "mutation", "effect", "FINAL CONFIDENCE GRADING"],
    )
    who = who[who["drug"] == drug].copy()
    who["truth"] = who["FINAL CONFIDENCE GRADING"].map(truth_label)
    if truth_mode == "resistant":
        who = who[who["truth"] == "R"].copy()
    elif truth_mode == "sensitive":
        who = who[who["truth"] == "S"].copy()
    who["mutation"] = who["mutation"].astype(str).str.strip()
    who["strict_comparable"] = who.apply(
        lambda row: is_strict_comparable_cds_variant(row["mutation"], row["effect"]),
        axis=1,
    )
    current_support = who.apply(
        lambda row: classify_current_support(
            row["gene"], row["mutation"], row["effect"], gene_window_specs
        ),
        axis=1,
    )
    who["current_supported"] = [bool(item[0]) for item in current_support]
    who["current_gap_reason"] = [item[1] for item in current_support]
    return who


def _simplify_db_record(row: sqlite3.Row) -> Dict[str, object]:
    return {
        "logic_type": str(row["logic_type"] or "POSITIVE").upper(),
        "status": str(row["status"] or ""),
        "final_interpretation": str(row["final_interpretation"] or ""),
        "resistance_phenotype": str(row["resistance_phenotype"] or ""),
    }


def load_db_status_by_key(db_path: Path, modeled_genes: Iterable[str]) -> Tuple[
    Dict[Tuple[str, str], List[Dict[str, object]]],
    Dict[Tuple[str, str, int, str, str], List[Dict[str, object]]],
]:
    modeled_genes = list(modeled_genes)
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    try:
        rows = conn.execute(
            "\n            SELECT g.name AS gene, g.logic_type, m.pos, m.ref, m.alt, m.aa_change,\n                     m.region_type, m.region_name, m.status, m.final_interpretation,\n                     m.resistance_phenotype\n            FROM mutations m\n            JOIN genes g ON g.id = m.gene_id\n            WHERE g.name IN ({})\n            ".format(
                ",".join(("?" for _ in modeled_genes))
            ),
            modeled_genes,
        ).fetchall()
    finally:
        conn.close()
    status_by_key: Dict[Tuple[str, str], List[Dict[str, object]]] = defaultdict(list)
    regulatory_status_by_key: Dict[
        Tuple[str, str, int, str, str], List[Dict[str, object]]
    ] = defaultdict(list)
    for row in rows:
        gene = str(row["gene"])
        region_type = str(row["region_type"] or "CDS").upper()
        record = _simplify_db_record(row)
        if region_type == "REGULATORY":
            regulatory_status_by_key[
                gene,
                str(row["region_name"] or "REGULATORY"),
                int(row["pos"]),
                str(row["ref"]),
                str(row["alt"]),
            ].append(record)
            continue
        status_by_key[gene, f"c.{row['pos']}{row['ref']}>{row['alt']}"].append(record)
        for mutation in aa_change_to_mutations(row["aa_change"]):
            status_by_key[gene, mutation].append(record)
    return (status_by_key, regulatory_status_by_key)


def derive_prediction_label(records: Iterable[Dict[str, object]]) -> Optional[str]:
    labels: Set[str] = set()
    for record in records:
        status = str(record.get("status") or "")
        logic_type = str(record.get("logic_type") or "POSITIVE").upper()
        interpretation = str(record.get("final_interpretation") or "").strip().lower()
        stored_phenotype = str(record.get("resistance_phenotype") or "").strip().upper()
        if stored_phenotype == "RESISTANT":
            labels.add("R")
            continue
        if stored_phenotype == "SENSITIVE":
            labels.add("S")
            continue
        if stored_phenotype == "UNKNOWN":
            pass
        if status == "RE_SENSITIZED":
            labels.add("S")
            continue
        if interpretation.startswith("sensitive:"):
            labels.add("S")
            continue
        if status == "LOSS_OF_FUNCTION":
            labels.add("S" if logic_type == "STRUCTURAL" else "R")
            continue
        if status == "COMPLETED":
            if logic_type in {"NEGATIVE", "STRUCTURAL"}:
                labels.add("S")
                continue
            if interpretation.startswith("resistant:"):
                labels.add("R")
                continue
    if len(labels) == 1:
        return next(iter(labels))
    return None


def summarize_scope(
    who: pd.DataFrame,
    status_by_key: Dict[Tuple[str, str], List[Dict[str, object]]],
    regulatory_status_by_key: Dict[
        Tuple[str, str, int, str, str], List[Dict[str, object]]
    ],
    gene_window_specs: Dict[str, List[Dict[str, Any]]],
    scope_column: str,
    modeled_genes: Iterable[str],
) -> Dict[str, object]:
    modeled_genes = set(modeled_genes)
    who_modeled = who[who["gene"].isin(modeled_genes)].copy()
    comparable = who_modeled[
        who_modeled["truth"].notna() & who_modeled[scope_column]
    ].copy()
    coverage = Counter()
    matched_final: List[Dict[str, object]] = []
    per_gene = defaultdict(Counter)
    for _, row in comparable.iterrows():
        gene = str(row["gene"])
        mutation = str(row["mutation"])
        truth = str(row["truth"])
        effect = str(row["effect"] or "")
        if effect == "upstream_gene_variant":
            records: List[Dict[str, object]] = []
            for key in regulatory_candidate_keys(gene, mutation, gene_window_specs):
                records.extend(regulatory_status_by_key.get(key, []))
        else:
            records = list(status_by_key.get((gene, mutation), []))
        statuses = Counter((str(record.get("status") or "") for record in records))
        pred = derive_prediction_label(records)
        per_gene[gene]["who_scope"] += 1
        if pred in {"R", "S"}:
            coverage["final"] += 1
            per_gene[gene]["final"] += 1
            matched_final.append(
                {
                    "gene": gene,
                    "mutation": mutation,
                    "truth": truth,
                    "pred": pred,
                    "grading": row["FINAL CONFIDENCE GRADING"],
                    "statuses": dict(statuses),
                    "records": records,
                }
            )
            if pred == truth:
                per_gene[gene]["correct"] += 1
            else:
                per_gene[gene]["wrong"] += 1
        elif "BOLTZ_READY" in statuses:
            coverage["pending_boltz"] += 1
            per_gene[gene]["pending_boltz"] += 1
        elif statuses:
            coverage["other"] += 1
            per_gene[gene]["other"] += 1
        else:
            coverage["absent"] += 1
            per_gene[gene]["absent"] += 1
    cm = Counter(((item["truth"], item["pred"]) for item in matched_final))
    tp = cm["R", "R"]
    tn = cm["S", "S"]
    fp = cm["S", "R"]
    fn = cm["R", "S"]
    total = tp + tn + fp + fn
    return {
        "who_total_drug": int(len(who)),
        "who_total_modeled": int(len(who_modeled)),
        "who_scope": int(len(comparable)),
        "who_gene_counts": who_modeled["gene"].value_counts().to_dict(),
        "coverage": dict(coverage),
        "matched_final": matched_final,
        "confusion": {"TP": tp, "TN": tn, "FP": fp, "FN": fn},
        "metrics": {
            "accuracy": round((tp + tn) / total, 4) if total else None,
            "precision": round(tp / (tp + fp), 4) if tp + fp else None,
            "recall": round(tp / (tp + fn), 4) if tp + fn else None,
            "specificity": round(tn / (tn + fp), 4) if tn + fp else None,
        },
        "per_gene": {
            gene: {
                "who_scope": counts.get("who_scope", 0),
                "final": counts.get("final", 0),
                "correct": counts.get("correct", 0),
                "wrong": counts.get("wrong", 0),
                "pending_boltz": counts.get("pending_boltz", 0),
                "other": counts.get("other", 0),
                "absent": counts.get("absent", 0),
                "final_accuracy": (
                    round(counts.get("correct", 0) / counts.get("final", 1), 4)
                    if counts.get("final", 0)
                    else None
                ),
            }
            for gene, counts in sorted(per_gene.items())
        },
    }


def main() -> int:
    args = parse_args()
    modeled_genes = resolve_modeled_genes(args.extra_modeled_gene)
    gene_window_specs = normalize_window_specs(args.config)
    who = load_who_catalog(args.xlsx, args.drug, args.truth, gene_window_specs)
    status_by_key, regulatory_status_by_key = load_db_status_by_key(
        args.db, modeled_genes
    )
    strict_summary = summarize_scope(
        who,
        status_by_key,
        regulatory_status_by_key,
        gene_window_specs,
        "strict_comparable",
        modeled_genes,
    )
    summary = summarize_scope(
        who,
        status_by_key,
        regulatory_status_by_key,
        gene_window_specs,
        "current_supported",
        modeled_genes,
    )
    who_genes = set(who["gene"].dropna())
    modeled_genes = set(modeled_genes)
    who_not_modeled_counts = (
        who.loc[~who["gene"].isin(modeled_genes), "gene"].value_counts().to_dict()
    )
    current_gap_counts = Counter(
        (
            str(reason)
            for reason in who.loc[
                who["gene"].isin(modeled_genes)
                & who["truth"].notna()
                & ~who["current_supported"],
                "current_gap_reason",
            ]
            if reason
        )
    )
    print(f"drug={args.drug}")
    print(f"truth_filter={args.truth}")
    print(f"who_xlsx={args.xlsx}")
    print(f"db={args.db}")
    print(f"config={args.config}")
    print(f"who_total_drug={summary['who_total_drug']}")
    print(f"who_total_modeled={summary['who_total_modeled']}")
    print(f"who_comparable={summary['who_scope']}")
    print(f"who_strict_comparable={strict_summary['who_scope']}")
    print(f"who_current_support_gaps={dict(current_gap_counts)}")
    print(f"who_gene_counts={summary['who_gene_counts']}")
    print(f"modeled_not_in_who={sorted(modeled_genes - who_genes)}")
    print(f"who_not_modeled={sorted(who_genes - modeled_genes)}")
    print(f"who_not_modeled_counts={who_not_modeled_counts}")
    print(f"coverage={summary['coverage']}")
    print(f"strict_coverage={strict_summary['coverage']}")
    print(f"confusion={summary['confusion']}")
    print(f"strict_confusion={strict_summary['confusion']}")
    print(f"metrics={summary['metrics']}")
    print(f"strict_metrics={strict_summary['metrics']}")
    print("per_gene=")
    for gene, gene_summary in summary["per_gene"].items():
        print(f"  {gene}: {gene_summary}")
    print("strict_per_gene=")
    for gene, gene_summary in strict_summary["per_gene"].items():
        print(f"  {gene}: {gene_summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
