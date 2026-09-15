from __future__ import annotations
import argparse
import csv
import re
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path
from typing import (
    TYPE_CHECKING,
    Any,
    Dict,
    Iterable,
    List,
    Mapping,
    Optional,
    Sequence,
    Set,
    Tuple,
)

if TYPE_CHECKING:
    import pandas as pd
else:
    try:
        import pandas as pd
    except ImportError:
        pd = None
try:
    import yaml
except ImportError:
    yaml = None
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
DEFAULT_PRIMARY_GENES = (
    "atpE",
    "Rv0678",
    "mmpL5",
    "mmpS5",
    "pepQ",
    "mtrA",
    "mtrB",
    "lpqB",
)
DEFAULT_ACCOUNTING_ONLY_GENES = ()
DEFAULT_EXPLORATORY_GENES = ("glpK", "Rv1979c")
SIMPLE_CDS_RE = re.compile("c\\.(\\d+)([ACGT]+)>([ACGT]+)$")
UPSTREAM_RE = re.compile("c\\.-(\\d+)([ACGT]+)>([ACGT]+)$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export publication-ready per-row WHO BDQ benchmark tables."
    )
    parser.add_argument(
        "--who-xlsx", type=Path, required=True, help="WHO catalogue workbook path."
    )
    parser.add_argument(
        "--db", type=Path, required=True, help="AMR-Hunter SQLite database path."
    )
    parser.add_argument(
        "--drug", default="Bedaquiline", help="Drug name in the WHO workbook."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("reports/publication"),
        help="Directory for publication CSV outputs.",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="Optional AMR-Hunter config.yaml used to resolve regulatory windows.",
    )
    parser.add_argument(
        "--primary-gene",
        action="append",
        default=[],
        help="Override/add a primary benchmark gene. Repeatable.",
    )
    parser.add_argument(
        "--accounting-only-gene",
        action="append",
        default=[],
        help="Override/add an accounting-only gene. Repeatable.",
    )
    return parser.parse_args()


def grade_number(grading: object) -> Optional[int]:
    text = str(grading or "").strip()
    match = re.match("(\\d)\\)", text)
    if match:
        return int(match.group(1))
    if "Assoc. w. R - interim" in text:
        return 2
    if "Assoc. w. R" in text:
        return 1
    if "Assoc. w. S - interim" in text:
        return 4
    if "Assoc. w. S" in text:
        return 5
    return None


def truth_label(grading: object) -> str:
    grade = grade_number(grading)
    if grade in {1, 2}:
        return "R"
    if grade in {4, 5}:
        return "S"
    return ""


def mutation_class(mutation: object, effect: object) -> str:
    mut = str(mutation or "").strip()
    eff = str(effect or "").strip()
    lower = mut.lower()
    if eff == "upstream_gene_variant" or UPSTREAM_RE.fullmatch(mut):
        return "regulatory"
    if mut == "LoF" or lower == "lof":
        return "generic_lof"
    if "fs" in lower or "frameshift" in lower:
        return "frameshift"
    if "del" in lower or "ins" in lower or "dup" in lower:
        return "indel"
    if mut.startswith("p.Met1") or mut.endswith("?"):
        return "start_lost"
    if mut.endswith("*") or mut.endswith("Ter"):
        return "nonsense"
    if SIMPLE_CDS_RE.fullmatch(mut):
        return "snv"
    if mut.startswith("p."):
        return "protein_substitution"
    return "unsupported"


def benchmark_role(
    gene: str,
    primary_genes: Set[str],
    accounting_only_genes: Set[str],
    exploratory_genes: Set[str],
) -> str:
    if gene in primary_genes:
        return "primary"
    if gene in accounting_only_genes:
        return "primary_accounting_only"
    if gene in exploratory_genes:
        return "exploratory"
    return "out_of_scope"


def normalize_window_specs(
    config_path: Optional[Path],
) -> Dict[str, List[Dict[str, Any]]]:
    if config_path is None or not config_path.exists():
        return {}
    if yaml is None:
        raise RuntimeError("PyYAML is required when --config is provided.")
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


def aa_change_to_mutations(change: object) -> Set[str]:
    text = str(change or "").strip()
    if not text:
        return set()
    if text.startswith("p."):
        if text.endswith("Ter"):
            return {text, f"{text[:-3]}*"}
        if text.endswith("*"):
            return {text, f"{text[:-1]}Ter"}
        return {text}
    if len(text) < 3:
        return set()
    wt, mt, pos = (text[0], text[-1], text[1:-1])
    if wt not in AA3 or mt not in AA3 or (not pos.isdigit()):
        return set()
    mutation = f"p.{AA3[wt]}{int(pos)}{AA3[mt]}"
    if mt == "*":
        return {mutation, f"p.{AA3[wt]}{int(pos)}Ter"}
    return {mutation}


def load_who_rows(xlsx_path: Path, drug: str) -> List[Dict[str, str]]:
    if pd is None:
        raise RuntimeError("pandas and openpyxl are required to read the WHO workbook.")
    who = pd.read_excel(
        xlsx_path,
        sheet_name="Catalogue_master_file",
        header=2,
        usecols=["drug", "gene", "mutation", "effect", "FINAL CONFIDENCE GRADING"],
    )
    who = who[who["drug"] == drug].copy()
    rows: List[Dict[str, str]] = []
    for idx, row in who.reset_index(drop=True).iterrows():
        grade = grade_number(row["FINAL CONFIDENCE GRADING"])
        rows.append(
            {
                "who_row_id": str(idx + 1),
                "drug": str(row["drug"] or ""),
                "gene": str(row["gene"] or "").strip(),
                "mutation": str(row["mutation"] or "").strip(),
                "effect": str(row["effect"] or "").strip(),
                "who_grade": f"G{grade}" if grade else "",
                "who_grading": str(row["FINAL CONFIDENCE GRADING"] or ""),
                "who_truth": truth_label(row["FINAL CONFIDENCE GRADING"]),
                "mutation_class": mutation_class(row["mutation"], row["effect"]),
            }
        )
    return rows


def simplify_record(row: sqlite3.Row) -> Dict[str, str]:
    return {
        "gene": str(row["gene"] or ""),
        "logic_type": str(row["logic_type"] or ""),
        "mutation": str(row["mutation_key"] or ""),
        "aa_change": str(row["aa_change"] or ""),
        "region_type": str(row["region_type"] or "CDS"),
        "region_name": str(row["region_name"] or ""),
        "status": str(row["status"] or ""),
        "resistance_phenotype": str(row["resistance_phenotype"] or ""),
        "functional_state": str(row["functional_state"] or ""),
        "evidence_code": str(row["evidence_code"] or ""),
        "interpretation": str(row["final_interpretation"] or ""),
    }


def load_db_records(
    db_path: Path,
) -> Tuple[
    Dict[Tuple[str, str], List[Dict[str, str]]],
    Dict[Tuple[str, str, int, str, str], List[Dict[str, str]]],
    Dict[Tuple[str, int, str, str], List[Dict[str, str]]],
]:
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    try:
        columns = {
            row["name"]
            for row in conn.execute("PRAGMA table_info(mutations)").fetchall()
        }
        phenotype_expr = (
            "m.resistance_phenotype"
            if "resistance_phenotype" in columns
            else "NULL AS resistance_phenotype"
        )
        functional_expr = (
            "m.functional_state"
            if "functional_state" in columns
            else "NULL AS functional_state"
        )
        evidence_expr = (
            "m.evidence_code" if "evidence_code" in columns else "NULL AS evidence_code"
        )
        rows = conn.execute(
            f"\n            SELECT\n                g.name AS gene,\n                g.logic_type AS logic_type,\n                m.pos,\n                m.ref,\n                m.alt,\n                'c.' || m.pos || m.ref || '>' || m.alt AS mutation_key,\n                m.aa_change,\n                m.region_type,\n                m.region_name,\n                m.status,\n                m.final_interpretation,\n                {phenotype_expr},\n                {functional_expr},\n                {evidence_expr}\n            FROM mutations m\n            JOIN genes g ON g.id = m.gene_id\n            "
        ).fetchall()
    finally:
        conn.close()
    cds: Dict[Tuple[str, str], List[Dict[str, str]]] = defaultdict(list)
    regulatory: Dict[Tuple[str, str, int, str, str], List[Dict[str, str]]] = (
        defaultdict(list)
    )
    regulatory_direct: Dict[Tuple[str, int, str, str], List[Dict[str, str]]] = (
        defaultdict(list)
    )
    for row in rows:
        record = simplify_record(row)
        gene = record["gene"]
        region_type = record["region_type"].upper()
        if region_type == "REGULATORY":
            region_name = str(row["region_name"] or "REGULATORY")
            pos = int(row["pos"])
            ref = str(row["ref"])
            alt = str(row["alt"])
            regulatory[gene, region_name, pos, ref, alt].append(record)
            regulatory_direct[gene, pos, ref, alt].append(record)
            continue
        cds[gene, record["mutation"]].append(record)
        for mut in aa_change_to_mutations(record["aa_change"]):
            cds[gene, mut].append(record)
    return (cds, regulatory, regulatory_direct)


def derive_prediction(records: Sequence[Mapping[str, str]]) -> Tuple[str, str]:
    labels: Set[str] = set()
    evidence_codes: List[str] = []
    for record in records:
        evidence_code = str(record.get("evidence_code") or "")
        if evidence_code:
            evidence_codes.append(evidence_code)
        phenotype = str(record.get("resistance_phenotype") or "").strip().upper()
        status = str(record.get("status") or "").strip().upper()
        logic_type = str(record.get("logic_type") or "POSITIVE").strip().upper()
        interpretation = str(record.get("interpretation") or "").strip().lower()
        if phenotype == "RESISTANT":
            labels.add("R")
            continue
        if phenotype == "SENSITIVE":
            labels.add("S")
            continue
        if status == "RE_SENSITIZED":
            labels.add("S")
            continue
        if interpretation.startswith("sensitive:"):
            labels.add("S")
            continue
        if interpretation.startswith("resistant:"):
            labels.add("R")
            continue
        if status == "LOSS_OF_FUNCTION":
            labels.add("S" if logic_type == "STRUCTURAL" else "R")
            continue
        if status == "COMPLETED" and logic_type in {"NEGATIVE", "STRUCTURAL"}:
            labels.add("S")
            continue
    evidence_summary = "|".join(sorted(set(evidence_codes)))
    if len(labels) == 1:
        return (next(iter(labels)), evidence_summary)
    if len(labels) > 1:
        return ("CONFLICT", evidence_summary)
    return ("", evidence_summary)


def regulatory_candidate_keys(
    gene: str, mutation: str, window_specs: Mapping[str, Sequence[Mapping[str, Any]]]
) -> Set[Tuple[str, str, int, str, str]]:
    match = UPSTREAM_RE.fullmatch(mutation)
    if not match:
        return set()
    upstream_distance, ref, alt = match.groups()
    distance = int(upstream_distance)
    candidate_keys: Set[Tuple[str, str, int, str, str]] = set()
    for spec in window_specs.get(gene, []):
        mode = str(spec.get("mode") or "")
        region_name = str(spec.get("name") or "REGULATORY")
        if mode == "upstream":
            length_bp = int(spec.get("length_bp") or 0)
            if 1 <= distance <= length_bp:
                candidate_keys.add(
                    (gene, region_name, length_bp - distance + 1, ref, alt)
                )
        elif mode == "around_start":
            upstream_bp = int(spec.get("upstream_bp") or 0)
            if 1 <= distance <= upstream_bp:
                candidate_keys.add(
                    (gene, region_name, upstream_bp - distance + 1, ref, alt)
                )
    return candidate_keys


def regulatory_direct_key(mutation: str) -> Optional[Tuple[int, str, str]]:
    match = UPSTREAM_RE.fullmatch(mutation)
    if not match:
        return None
    position, ref, alt = match.groups()
    return (int(position), ref, alt)


def find_records(
    who_row: Mapping[str, str],
    cds: Mapping[Tuple[str, str], List[Dict[str, str]]],
    regulatory: Mapping[Tuple[str, str, int, str, str], List[Dict[str, str]]],
    regulatory_direct: Mapping[Tuple[str, int, str, str], List[Dict[str, str]]],
    window_specs: Mapping[str, Sequence[Mapping[str, Any]]],
) -> List[Dict[str, str]]:
    gene = who_row["gene"]
    mutation = who_row["mutation"]
    records = list(cds.get((gene, mutation), []))
    candidate_keys = regulatory_candidate_keys(gene, mutation, window_specs)
    for key in candidate_keys:
        records.extend(regulatory.get(key, []))
    if not candidate_keys and gene not in window_specs:
        direct_key = regulatory_direct_key(mutation)
        if direct_key is not None:
            pos, ref, alt = direct_key
            records.extend(regulatory_direct.get((gene, pos, ref, alt), []))
    return records


def current_evaluator_scope_gap(
    row: Mapping[str, str],
    role: str,
    window_specs: Mapping[str, Sequence[Mapping[str, Any]]],
) -> str:
    if not row["who_truth"] or role != "primary":
        return ""
    if row["mutation_class"] == "generic_lof":
        return "generic_lof_category"
    if (
        row["mutation_class"] == "regulatory"
        and row["effect"] == "upstream_gene_variant"
    ):
        if window_specs and (
            not regulatory_candidate_keys(row["gene"], row["mutation"], window_specs)
        ):
            return "regulatory_out_of_config_window"
    return ""


def gap_reason(
    row: Mapping[str, str],
    role: str,
    records: Sequence[Mapping[str, str]],
    pred: str,
    scope_gap: str,
) -> str:
    if not row["who_truth"]:
        return "excluded_uncertain_g3_or_unmapped_grade"
    if role == "out_of_scope":
        return "gene_out_of_benchmark_scope"
    if role == "exploratory":
        return "exploratory_gene_excluded_from_primary_denominator"
    if scope_gap:
        return scope_gap
    if pred in {"R", "S"}:
        return ""
    if pred == "CONFLICT":
        return "conflicting_prediction_records"
    if records:
        return "records_without_final_rs_prediction"
    mut_class = row["mutation_class"]
    if mut_class == "generic_lof":
        return "generic_lof_category"
    if mut_class in {"frameshift", "indel", "start_lost"}:
        return "unsupported_non_snp_variant"
    if mut_class == "regulatory":
        return "regulatory_variant_not_matched"
    return "missing_amr_hunter_record"


def source_summary(records: Sequence[Mapping[str, str]]) -> str:
    parts = []
    for record in records[:3]:
        parts.append(
            ":".join(
                [
                    str(record.get("logic_type") or ""),
                    str(record.get("status") or ""),
                    str(record.get("resistance_phenotype") or ""),
                    str(record.get("evidence_code") or ""),
                ]
            )
        )
    if len(records) > 3:
        parts.append(f"+{len(records) - 3}_more")
    return ";".join(parts)


def status_summary(records: Sequence[Mapping[str, str]]) -> str:
    statuses = sorted(
        {str(record.get("status") or "") for record in records if record.get("status")}
    )
    return "|".join(statuses)


def build_rows(
    who_rows: Sequence[Dict[str, str]],
    cds: Mapping[Tuple[str, str], List[Dict[str, str]]],
    regulatory: Mapping[Tuple[str, str, int, str, str], List[Dict[str, str]]],
    regulatory_direct: Mapping[Tuple[str, int, str, str], List[Dict[str, str]]],
    window_specs: Mapping[str, Sequence[Mapping[str, Any]]],
    primary_genes: Set[str],
    accounting_only_genes: Set[str],
    exploratory_genes: Set[str],
) -> List[Dict[str, str]]:
    output = []
    for row in who_rows:
        role = benchmark_role(
            row["gene"], primary_genes, accounting_only_genes, exploratory_genes
        )
        records = find_records(row, cds, regulatory, regulatory_direct, window_specs)
        pred, evidence_codes = derive_prediction(records)
        scope_gap = current_evaluator_scope_gap(row, role, window_specs)
        gap = gap_reason(row, role, records, pred, scope_gap)
        correct = ""
        matched_status = "excluded"
        if row["who_truth"] and role == "primary":
            if pred in {"R", "S"}:
                matched_status = "matched_final"
                correct = "true" if pred == row["who_truth"] else "false"
            elif pred == "CONFLICT":
                matched_status = "conflict"
                correct = "false"
            else:
                matched_status = "unmatched"
                correct = "false"
        elif row["who_truth"] and role == "primary_accounting_only":
            matched_status = "accounting_only"
        elif row["who_truth"] and role == "exploratory":
            matched_status = "exploratory_excluded"
        output.append(
            {
                **row,
                "current_evaluator_scope": (
                    ("false" if scope_gap else "true")
                    if row["who_truth"] and role == "primary"
                    else ""
                ),
                "benchmark_role": role,
                "amr_hunter_prediction": pred,
                "matched_status": matched_status,
                "correct": correct,
                "gap_reason": gap,
                "record_count": str(len(records)),
                "source_statuses": status_summary(records),
                "evidence_code": evidence_codes,
                "source_records": source_summary(records),
            }
        )
    return output


def summarize(rows: Sequence[Mapping[str, str]]) -> List[Dict[str, str]]:
    buckets: Dict[Tuple[str, str], List[Mapping[str, str]]] = defaultdict(list)
    for row in rows:
        buckets["overall_primary", "all"].append(row)
        buckets["gene", row["gene"]].append(row)
        buckets["role", row["benchmark_role"]].append(row)
    summaries: List[Dict[str, str]] = []
    for (scope_type, scope_value), items in sorted(buckets.items()):
        evaluable = [
            item
            for item in items
            if item["who_truth"] in {"R", "S"} and item["benchmark_role"] == "primary"
        ]
        current_scope = [
            item for item in evaluable if item["current_evaluator_scope"] == "true"
        ]
        matched = [
            item for item in evaluable if item["matched_status"] == "matched_final"
        ]
        correct = [item for item in matched if item["correct"] == "true"]
        wrong = [item for item in matched if item["correct"] == "false"]
        unmatched = [
            item for item in evaluable if item["matched_status"] == "unmatched"
        ]
        cm = Counter(
            ((item["who_truth"], item["amr_hunter_prediction"]) for item in matched)
        )
        total_matched = len(matched)
        total_evaluable = len(evaluable)
        summaries.append(
            {
                "scope_type": scope_type,
                "scope_value": scope_value,
                "who_rows": str(len(items)),
                "evaluable_primary_rows": str(total_evaluable),
                "current_evaluator_scope_rows": str(len(current_scope)),
                "matched_rows": str(total_matched),
                "correct_rows": str(len(correct)),
                "wrong_rows": str(len(wrong)),
                "unmatched_rows": str(len(unmatched)),
                "coverage": (
                    f"{total_matched / total_evaluable:.4f}" if total_evaluable else ""
                ),
                "current_scope_coverage": (
                    f"{total_matched / len(current_scope):.4f}" if current_scope else ""
                ),
                "matched_call_accuracy": (
                    f"{len(correct) / total_matched:.4f}" if total_matched else ""
                ),
                "tp": str(cm["R", "R"]),
                "tn": str(cm["S", "S"]),
                "fp": str(cm["S", "R"]),
                "fn": str(cm["R", "S"]),
            }
        )
    return summaries


def write_csv(
    path: Path, rows: Sequence[Mapping[str, str]], fieldnames: Sequence[str]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def main() -> int:
    args = parse_args()
    primary_genes = set(DEFAULT_PRIMARY_GENES)
    accounting_only_genes = set(DEFAULT_ACCOUNTING_ONLY_GENES)
    exploratory_genes = set(DEFAULT_EXPLORATORY_GENES)
    if args.primary_gene:
        primary_genes.update(args.primary_gene)
    if args.accounting_only_gene:
        accounting_only_genes.update(args.accounting_only_gene)
    who_rows = load_who_rows(args.who_xlsx, args.drug)
    window_specs = normalize_window_specs(args.config)
    cds, regulatory, regulatory_direct = load_db_records(args.db)
    rows = build_rows(
        who_rows,
        cds,
        regulatory,
        regulatory_direct,
        window_specs,
        primary_genes,
        accounting_only_genes,
        exploratory_genes,
    )
    summaries = summarize(rows)
    per_row_fields = [
        "who_row_id",
        "drug",
        "gene",
        "mutation",
        "effect",
        "who_grade",
        "who_truth",
        "who_grading",
        "mutation_class",
        "current_evaluator_scope",
        "benchmark_role",
        "amr_hunter_prediction",
        "matched_status",
        "correct",
        "gap_reason",
        "record_count",
        "source_statuses",
        "evidence_code",
        "source_records",
    ]
    summary_fields = [
        "scope_type",
        "scope_value",
        "who_rows",
        "evaluable_primary_rows",
        "current_evaluator_scope_rows",
        "matched_rows",
        "correct_rows",
        "wrong_rows",
        "unmatched_rows",
        "coverage",
        "current_scope_coverage",
        "matched_call_accuracy",
        "tp",
        "tn",
        "fp",
        "fn",
    ]
    write_csv(args.output_dir / "who_bdq_per_row_benchmark.csv", rows, per_row_fields)
    write_csv(
        args.output_dir / "who_bdq_benchmark_summary.csv", summaries, summary_fields
    )
    print(f"per_row={args.output_dir / 'who_bdq_per_row_benchmark.csv'}")
    print(f"summary={args.output_dir / 'who_bdq_benchmark_summary.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
