from __future__ import annotations
import argparse
import csv
import sqlite3
from pathlib import Path
from core.analyzer import resolve_analysis_decision


def _build_mutation_label(row: sqlite3.Row) -> str:
    aa_change = row["aa_change"]
    if aa_change:
        return str(aa_change)
    region_name = row["region_name"] or row["region_type"] or "NONCODING"
    return f"{region_name}:{row['ref']}{row['pos']}{row['alt']}"


def _has_boltz_core(row: sqlite3.Row) -> bool:
    return any(
        (row[key] is not None for key in ("boltz_ptm", "boltz_iptm", "boltz_plddt"))
    )


def _derive_evidence_tier(row: sqlite3.Row) -> tuple[str, str, str]:
    status = str(row["status"] or "")
    logic_type = str(row["logic_type"] or "")
    has_boltz = _has_boltz_core(row)
    if status == "RE_SENSITIZED":
        return ("HIGH", "PATHWAY_RE_SENSITIZED", "sensitive_or_resensitized")
    if status == "LOSS_OF_FUNCTION":
        label = (
            "likely_resistant"
            if logic_type == "NEGATIVE"
            else "likely_functional_collapse"
        )
        return ("HIGH", "EVO_SHORT_CIRCUIT", label)
    if status == "LETHAL_SKIP":
        return ("HIGH", "EVO_LETHAL", "unlikely_resistance_lethal")
    if status == "FAILED":
        return ("LOW", "FAILED", "needs_retry")
    if has_boltz:
        return ("HIGH", "EVO_PLUS_BOLTZ", "structure_reviewed")
    if status == "COMPLETED":
        return ("MEDIUM", "EVO_ONLY_COMPLETED", "evo_only_completed")
    if status == "BOLTZ_READY":
        return ("MEDIUM", "QUEUED_FOR_BOLTZ", "awaiting_boltz")
    return ("LOW", "OTHER", "manual_review")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Export tiered AMR mutation analysis table"
    )
    parser.add_argument("--db-path", required=True, help="Path to the SQLite database")
    parser.add_argument("--output", required=True, help="Path to the CSV output file")
    args = parser.parse_args()
    conn = sqlite3.connect(Path(args.db_path))
    conn.row_factory = sqlite3.Row
    rows = conn.execute(
        "\n        SELECT g.name AS gene_name, g.logic_type,\n               m.id, m.aa_change, m.region_type, m.region_name, m.pos, m.ref, m.alt,\n               m.evo_model, m.evo_delta, m.status, m.final_interpretation,\n             m.resistance_phenotype, m.functional_state, m.evidence_code,\n               m.boltz_ptm, m.boltz_iptm, m.boltz_plddt,\n               m.boltz_stability_score, m.boltz_clash,\n               m.boltz_binding_affinity, m.boltz_affinity_pred_value,\n               m.boltz_pair_energy, m.boltz_complex_energy,\n               m.synergy_score, m.synergy_tags\n        FROM mutations m\n        JOIN genes g ON m.gene_id = g.id\n        ORDER BY g.name, m.id\n        "
    ).fetchall()
    conn.close()
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as fp:
        writer = csv.writer(fp)
        writer.writerow(
            [
                "gene",
                "logic_type",
                "mutation",
                "region_type",
                "evo_model",
                "evo_delta",
                "status",
                "evidence_tier",
                "evidence_source",
                "final_call",
                "boltz_core_available",
                "ptm",
                "iptm",
                "plddt",
                "stability_score",
                "clash",
                "aux_binding_affinity",
                "aux_affinity_pred_value",
                "aux_pair_energy",
                "aux_complex_energy",
                "synergy_score",
                "synergy_tags",
                "resistance_label",
                "functional_state",
                "evidence_code",
                "interpretation",
            ]
        )
        for row in rows:
            evidence_tier, evidence_source, final_call = _derive_evidence_tier(row)
            decision = resolve_analysis_decision(dict(row))
            writer.writerow(
                [
                    row["gene_name"],
                    row["logic_type"],
                    _build_mutation_label(row),
                    row["region_type"],
                    row["evo_model"] or "evo2_7b",
                    row["evo_delta"],
                    row["status"],
                    evidence_tier,
                    evidence_source,
                    final_call,
                    int(_has_boltz_core(row)),
                    row["boltz_ptm"],
                    row["boltz_iptm"],
                    row["boltz_plddt"],
                    row["boltz_stability_score"],
                    row["boltz_clash"],
                    row["boltz_binding_affinity"],
                    row["boltz_affinity_pred_value"],
                    row["boltz_pair_energy"],
                    row["boltz_complex_energy"],
                    row["synergy_score"],
                    row["synergy_tags"],
                    decision.resistance_phenotype,
                    decision.functional_state,
                    decision.evidence_code,
                    row["final_interpretation"],
                ]
            )
    print(f"exported_rows={len(rows)}")
    print(f"output={output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
