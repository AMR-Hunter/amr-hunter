from __future__ import annotations
import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

PRIMARY_BDQ_GENES = (
    "atpE",
    "Rv0678",
    "mmpL5",
    "mmpS5",
    "pepQ",
    "mtrA",
    "mtrB",
    "Rv1979c",
)
EXPLORATORY_BDQ_GENES = ("glpK",)
ACCOUNTING_ONLY_GENES = ("lpqB",)
TARGET_GENES = PRIMARY_BDQ_GENES + EXPLORATORY_BDQ_GENES + ACCOUNTING_ONLY_GENES
RESISTANCE_EVIDENCE_CODES = {
    "ATP_SYNTHASE_COMPENSATION",
    "NEGATIVE_LOF_REVIEWED",
    "NEGATIVE_LOF_SIGNATURE",
    "MTRAB_DEREGULATION",
}
SENSITIVITY_EVIDENCE_CODES = {
    "EFFLUX_PUMP_COLLAPSE",
    "NEGATIVE_FUNCTION_RETAINED",
    "SYNONYMOUS_NO_PROTEIN_CHANGE",
}
REGULATORY_EVIDENCE_CODES = {
    "REGULATORY_SCENARIO_UNSUPPORTED",
    "EFFLUX_PUMP_COLLAPSE",
    "MTRAB_DEREGULATION",
}
AA_MUTATION_RE = re.compile("^([A-Z*])(\\d+)([A-Z*])$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export publication-ready BDQ virtual mutation summary tables."
    )
    parser.add_argument(
        "--snapshot",
        type=Path,
        default=Path("reports/amr_results_snapshot.csv"),
        help="Existing AMR-Hunter result snapshot CSV.",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("config/config.yaml"),
        help="Optional config.yaml used only to annotate regulatory support.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("reports/publication"),
        help="Directory for publication CSV outputs.",
    )
    parser.add_argument(
        "--candidate-limit-per-gene",
        type=int,
        default=25,
        help="Maximum candidate/example rows exported for each gene and role.",
    )
    return parser.parse_args()


def read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(
    path: Path, rows: Sequence[Mapping[str, object]], fieldnames: Sequence[str]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def benchmark_role(gene: str) -> str:
    if gene in PRIMARY_BDQ_GENES:
        return "primary_bdq"
    if gene in EXPLORATORY_BDQ_GENES:
        return "exploratory_bdq"
    if gene in ACCOUNTING_ONLY_GENES:
        return "accounting_only_missing"
    return "out_of_scope"


def mutation_kind(mutation: str) -> str:
    if ":" in mutation or mutation.startswith(
        ("DEFAULT_", "UPSTREAM_", "RV", "MTR", "ATP")
    ):
        return "regulatory"
    match = AA_MUTATION_RE.fullmatch(mutation)
    if match:
        ref, _, alt = match.groups()
        if ref == alt:
            return "synonymous_or_same"
        if alt == "*":
            return "truncating"
        return "coding_substitution"
    return "other"


def numeric(value: object) -> Optional[float]:
    text = str(value or "").strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def compact_counts(counter: Counter, limit: int = 6) -> str:
    parts = []
    for key, count in counter.most_common(limit):
        parts.append(f"{key}:{count}")
    return ";".join(parts)


def paper_use_for_gene(gene: str, rows: Sequence[Mapping[str, str]]) -> str:
    if gene in ACCOUNTING_ONLY_GENES:
        return "not present in current snapshot; report as limitation"
    if gene in EXPLORATORY_BDQ_GENES:
        return "exploratory only; do not include in primary WHO denominator"
    if not rows:
        return "missing from current snapshot"
    return "main BDQ virtual-mutation discovery scope"


def summarize(rows: Sequence[Dict[str, str]]) -> List[Dict[str, object]]:
    by_gene: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_gene[row["gene"]].append(row)
    summary_rows: List[Dict[str, object]] = []
    all_genes = sorted(
        set(by_gene)
        | set(PRIMARY_BDQ_GENES)
        | set(EXPLORATORY_BDQ_GENES)
        | set(ACCOUNTING_ONLY_GENES)
    )
    for gene in all_genes:
        gene_rows = by_gene.get(gene, [])
        kind_counts = Counter((mutation_kind(row["mutation"]) for row in gene_rows))
        label_counts = Counter(
            (row["resistance_label"] or "BLANK" for row in gene_rows)
        )
        status_counts = Counter((row["final_status"] or "BLANK" for row in gene_rows))
        evidence_counts = Counter(
            (row["evidence_code"] or "BLANK" for row in gene_rows)
        )
        summary_rows.append(
            {
                "gene": gene,
                "benchmark_role": benchmark_role(gene),
                "total_rows": len(gene_rows),
                "coding_substitution_rows": kind_counts["coding_substitution"],
                "truncating_rows": kind_counts["truncating"],
                "synonymous_or_same_rows": kind_counts["synonymous_or_same"],
                "regulatory_rows": kind_counts["regulatory"],
                "other_rows": kind_counts["other"],
                "resistant_rows": label_counts["RESISTANT"],
                "sensitive_rows": label_counts["SENSITIVE"],
                "unknown_rows": label_counts["UNKNOWN"],
                "completed_rows": status_counts["COMPLETED"],
                "loss_of_function_rows": status_counts["LOSS_OF_FUNCTION"],
                "top_evidence_codes": compact_counts(evidence_counts),
                "paper_use": paper_use_for_gene(gene, gene_rows),
            }
        )
    return summary_rows


def row_priority(row: Mapping[str, str]) -> Tuple[int, float, str]:
    evidence = row["evidence_code"]
    kind = mutation_kind(row["mutation"])
    evidence_rank = {
        "ATP_SYNTHASE_COMPENSATION": 0,
        "NEGATIVE_LOF_SIGNATURE": 1,
        "NEGATIVE_LOF_REVIEWED": 2,
        "MTRAB_DEREGULATION": 3,
        "NEGATIVE_FUNCTION_RETAINED": 4,
        "EFFLUX_PUMP_COLLAPSE": 5,
        "SYNONYMOUS_NO_PROTEIN_CHANGE": 6,
    }.get(evidence, 9)
    kind_penalty = {
        "truncating": 0,
        "coding_substitution": 1,
        "other": 2,
        "regulatory": 3,
    }.get(kind, 4)
    delta = numeric(row.get("evo_delta"))
    magnitude = abs(delta) if delta is not None else 0.0
    return (evidence_rank + kind_penalty, -magnitude, row["mutation"])


def candidate_role(row: Mapping[str, str]) -> str:
    label = row["resistance_label"]
    evidence = row["evidence_code"]
    kind = mutation_kind(row["mutation"])
    if (
        label == "RESISTANT"
        and kind == "regulatory"
        and (evidence in RESISTANCE_EVIDENCE_CODES)
    ):
        return "regulatory_resistance_candidate"
    if (
        label == "RESISTANT"
        and evidence in RESISTANCE_EVIDENCE_CODES
        and (kind != "synonymous_or_same")
    ):
        return "resistance_candidate"
    if (
        label == "SENSITIVE"
        and evidence in SENSITIVITY_EVIDENCE_CODES
        and (kind != "synonymous_or_same")
    ):
        return "sensitivity_or_retained_function_example"
    if (
        label == "UNKNOWN"
        and evidence == "REGULATORY_SCENARIO_UNSUPPORTED"
        and (kind == "regulatory")
    ):
        return "unsupported_regulatory_evo_only"
    return ""


def candidate_note(row: Mapping[str, str]) -> str:
    role = candidate_role(row)
    evidence = row["evidence_code"]
    if role == "resistance_candidate":
        if evidence == "ATP_SYNTHASE_COMPENSATION":
            return "ATP synthase candidate; interpret as follow-up priority, not a clinical call"
        if evidence in {"NEGATIVE_LOF_SIGNATURE", "NEGATIVE_LOF_REVIEWED"}:
            return (
                "negative-logic LoF candidate; supports mechanism-driven prioritization"
            )
        if evidence == "MTRAB_DEREGULATION":
            return "mtrAB pathway candidate; useful as secondary mechanism panel"
    if role == "regulatory_resistance_candidate":
        return "regulatory/intergenic resistance candidate; use as curated extensibility example, not broad regulatory validation"
    if role == "sensitivity_or_retained_function_example":
        return "mechanistic sensitive/control example for explaining rule direction"
    if role == "unsupported_regulatory_evo_only":
        return "regulatory Evo-only result; useful for limitations, not primary claims"
    return ""


def candidate_rows(
    rows: Sequence[Dict[str, str]], limit_per_gene: int
) -> List[Dict[str, object]]:
    grouped: Dict[Tuple[str, str], List[Dict[str, str]]] = defaultdict(list)
    for row in rows:
        role = candidate_role(row)
        if not role:
            continue
        if benchmark_role(row["gene"]) not in {"primary_bdq", "exploratory_bdq"}:
            continue
        grouped[row["gene"], role].append(row)
    output: List[Dict[str, object]] = []
    for (gene, role), items in sorted(grouped.items()):
        seen = set()
        kept = 0
        for row in sorted(items, key=row_priority):
            key = (
                row["gene"],
                row["mutation"],
                row["evidence_code"],
                row["resistance_label"],
            )
            if key in seen:
                continue
            seen.add(key)
            output.append(
                {
                    "gene": gene,
                    "benchmark_role": benchmark_role(gene),
                    "candidate_role": role,
                    "mutation": row["mutation"],
                    "mutation_kind": mutation_kind(row["mutation"]),
                    "resistance_label": row["resistance_label"],
                    "functional_state": row["functional_state"],
                    "evidence_code": row["evidence_code"],
                    "final_status": row["final_status"],
                    "evo_delta": row["evo_delta"],
                    "stability_score": row["stability_score"],
                    "ptm": row["ptm"],
                    "iptm": row["iptm"],
                    "plddt": row["plddt"],
                    "synergy_score": row["synergy_score"],
                    "synergy_tags": row["synergy_tags"],
                    "paper_note": candidate_note(row),
                    "interpretation_excerpt": str(row["interpretation"] or "")[:240],
                }
            )
            kept += 1
            if kept >= limit_per_gene:
                break
    return output


def parse_config_gene_scenarios(config_path: Path) -> Dict[str, Dict[str, str]]:
    if not config_path.exists():
        return {}
    meta: Dict[str, Dict[str, str]] = {}
    current_gene = ""
    in_target_genes = False
    for raw_line in config_path.read_text(
        encoding="utf-8", errors="replace"
    ).splitlines():
        line = raw_line.strip()
        if line == "target_genes:":
            in_target_genes = True
            continue
        if in_target_genes and raw_line and (not raw_line.startswith((" ", "-"))):
            break
        if not in_target_genes:
            continue
        if line.startswith("gene:"):
            current_gene = line.split(":", 1)[1].strip().strip('"')
            meta.setdefault(current_gene, {})
        elif current_gene and line.startswith("scenario:"):
            meta[current_gene]["scenario"] = line.split(":", 1)[1].strip().strip('"')
        elif current_gene and line.startswith("include_upstream:"):
            meta[current_gene]["include_upstream"] = (
                line.split(":", 1)[1].strip().strip('"')
            )
        elif (
            current_gene
            and line.startswith("upstream_bp:")
            and (meta[current_gene].get("include_upstream") == "true")
            and ("include_upstream_bp" not in meta[current_gene])
        ):
            meta[current_gene]["include_upstream_bp"] = (
                line.split(":", 1)[1].strip().strip('"')
            )
        elif current_gene and line.startswith("- name:"):
            window = line.split(":", 1)[1].strip().strip('"')
            existing = meta[current_gene].get("noncoding_windows", "")
            meta[current_gene]["noncoding_windows"] = ";".join(
                [item for item in [existing, window] if item]
            )
    return meta


def regulatory_assessment_rows(
    rows: Sequence[Dict[str, str]], config_path: Path
) -> List[Dict[str, object]]:
    config_meta = parse_config_gene_scenarios(config_path)
    by_gene: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for row in rows:
        if mutation_kind(row["mutation"]) == "regulatory":
            by_gene[row["gene"]].append(row)
    output: List[Dict[str, object]] = []
    genes = sorted(set(by_gene) | set(PRIMARY_BDQ_GENES) | set(EXPLORATORY_BDQ_GENES))
    for gene in genes:
        gene_rows = by_gene.get(gene, [])
        label_counts = Counter(
            (row["resistance_label"] or "BLANK" for row in gene_rows)
        )
        evidence_counts = Counter(
            (row["evidence_code"] or "BLANK" for row in gene_rows)
        )
        meta = config_meta.get(gene, {})
        if gene == "Rv0678" and (not gene_rows):
            recommendation = "configured for PROTEIN_DNA/operator logic but no regulatory rows in current snapshot; run a small targeted regulatory job only if the manuscript needs this case study"
        elif gene == "Rv0678":
            recommendation = "best curated regulatory/intergenic extensibility example; present as mechanism-supported candidate evidence, not broad regulatory validation"
        elif label_counts["SENSITIVE"] and "EFFLUX_PUMP_COLLAPSE" in evidence_counts:
            recommendation = (
                "usable conservative noncoding example for sensitivity/control behavior"
            )
        elif label_counts["RESISTANT"] and "MTRAB_DEREGULATION" in evidence_counts:
            recommendation = (
                "secondary regulatory/pathway example; not Rv0678 protein-DNA binding"
            )
        elif (
            label_counts["UNKNOWN"]
            and "REGULATORY_SCENARIO_UNSUPPORTED" in evidence_counts
        ):
            recommendation = (
                "use as limitation/evo-only regulatory scan, not as resistance claim"
            )
        else:
            recommendation = "not recommended as a main regulatory result"
        output.append(
            {
                "gene": gene,
                "benchmark_role": benchmark_role(gene),
                "configured_scenario": meta.get("scenario", ""),
                "include_upstream": meta.get("include_upstream", ""),
                "include_upstream_bp": meta.get("include_upstream_bp", ""),
                "noncoding_windows": meta.get("noncoding_windows", ""),
                "snapshot_regulatory_rows": len(gene_rows),
                "regulatory_resistant_rows": label_counts["RESISTANT"],
                "regulatory_sensitive_rows": label_counts["SENSITIVE"],
                "regulatory_unknown_rows": label_counts["UNKNOWN"],
                "regulatory_evidence_codes": compact_counts(evidence_counts),
                "paper_recommendation": recommendation,
            }
        )
    return output


def main() -> int:
    args = parse_args()
    rows = read_csv(args.snapshot)
    rows = [row for row in rows if row.get("gene") in TARGET_GENES]
    summary = summarize(rows)
    candidates = candidate_rows(rows, args.candidate_limit_per_gene)
    regulatory = regulatory_assessment_rows(rows, args.config)
    summary_fields = [
        "gene",
        "benchmark_role",
        "total_rows",
        "coding_substitution_rows",
        "truncating_rows",
        "synonymous_or_same_rows",
        "regulatory_rows",
        "other_rows",
        "resistant_rows",
        "sensitive_rows",
        "unknown_rows",
        "completed_rows",
        "loss_of_function_rows",
        "top_evidence_codes",
        "paper_use",
    ]
    candidate_fields = [
        "gene",
        "benchmark_role",
        "candidate_role",
        "mutation",
        "mutation_kind",
        "resistance_label",
        "functional_state",
        "evidence_code",
        "final_status",
        "evo_delta",
        "stability_score",
        "ptm",
        "iptm",
        "plddt",
        "synergy_score",
        "synergy_tags",
        "paper_note",
        "interpretation_excerpt",
    ]
    regulatory_fields = [
        "gene",
        "benchmark_role",
        "configured_scenario",
        "include_upstream",
        "include_upstream_bp",
        "noncoding_windows",
        "snapshot_regulatory_rows",
        "regulatory_resistant_rows",
        "regulatory_sensitive_rows",
        "regulatory_unknown_rows",
        "regulatory_evidence_codes",
        "paper_recommendation",
    ]
    write_csv(args.output_dir / "virtual_mutation_summary.csv", summary, summary_fields)
    write_csv(
        args.output_dir / "virtual_mutation_candidates.csv",
        candidates,
        candidate_fields,
    )
    write_csv(
        args.output_dir / "regulatory_example_assessment.csv",
        regulatory,
        regulatory_fields,
    )
    print(f"summary={args.output_dir / 'virtual_mutation_summary.csv'}")
    print(f"candidates={args.output_dir / 'virtual_mutation_candidates.csv'}")
    print(f"regulatory={args.output_dir / 'regulatory_example_assessment.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
