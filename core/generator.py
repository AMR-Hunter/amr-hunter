from __future__ import annotations
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional
import yaml
from Bio.Seq import Seq
from core.database import DatabaseManager


class MutationEngine:

    def __init__(
        self,
        config_path: str = "config/config.yaml",
        db_manager: Optional[DatabaseManager] = None,
        codon_table: int = 11,
    ) -> None:
        self.config_path = Path(config_path)
        self.config = self._load_config()
        self.db = db_manager
        self.codon_table = codon_table

    def _load_config(self) -> Dict[str, Any]:
        with self.config_path.open("r", encoding="utf-8") as file:
            return yaml.safe_load(file) or {}

    def get_reference_genome_path(self, species_name: str) -> Optional[str]:
        for item in self.config.get("reference_genomes", []):
            if item.get("species") == species_name:
                return item.get("path")
        return None

    def generate_saturation_mutagenesis(
        self, gene_record: Dict[str, Any]
    ) -> List[Dict[str, Any]]:
        gene_id = int(gene_record["gene_id"])
        sequence = str(gene_record["sequence"]).upper().replace("U", "T")
        target_drug_id = gene_record.get("target_drug_id")
        return self._generate_region_mutations(
            gene_id=gene_id,
            sequence=sequence,
            target_drug_id=target_drug_id,
            region_type="CDS",
            region_name=None,
            codon_aware=True,
        )

    def generate_regulatory_mutagenesis(
        self, gene_record: Dict[str, Any], region_name: str = "upstream"
    ) -> List[Dict[str, Any]]:
        gene_id = int(gene_record["gene_id"])
        sequence = str(gene_record["sequence"]).upper().replace("U", "T")
        target_drug_id = gene_record.get("target_drug_id")
        return self._generate_region_mutations(
            gene_id=gene_id,
            sequence=sequence,
            target_drug_id=target_drug_id,
            region_type="REGULATORY",
            region_name=region_name,
            codon_aware=False,
        )

    def _generate_region_mutations(
        self,
        gene_id: int,
        sequence: str,
        target_drug_id: Optional[int],
        region_type: str,
        region_name: Optional[str],
        codon_aware: bool,
    ) -> List[Dict[str, Any]]:
        if not sequence:
            return []
        if codon_aware and len(sequence) % 3 != 0:
            raise ValueError(f"Gene {gene_id} CDS length is not multiple of 3.")
        nucleotides = ("A", "C", "G", "T")
        wt_aa_seq = (
            str(Seq(sequence).translate(table=self.codon_table, to_stop=False))
            if codon_aware
            else ""
        )
        mutations: List[Dict[str, Any]] = []
        for idx, ref_nt in enumerate(sequence):
            if ref_nt not in nucleotides:
                continue
            codon_start = idx // 3 * 3
            codon_offset = idx % 3
            wt_codon = sequence[codon_start : codon_start + 3]
            if any((base not in nucleotides for base in wt_codon)):
                continue
            wt_aa = wt_aa_seq[idx // 3] if codon_aware else ""
            aa_pos = idx // 3 + 1 if codon_aware else 0
            for alt_nt in nucleotides:
                if alt_nt == ref_nt:
                    continue
                mt_codon = list(wt_codon)
                mt_codon[codon_offset] = alt_nt
                mt_codon_str = "".join(mt_codon)
                aa_change = None
                mutation_type = "NONCODING"
                if codon_aware:
                    mt_aa = str(
                        Seq(mt_codon_str).translate(
                            table=self.codon_table, to_stop=False
                        )
                    )
                    aa_change = f"{wt_aa}{aa_pos}{mt_aa}"
                    mutation_type = "SYNONYMOUS" if wt_aa == mt_aa else "NONSYNONYMOUS"
                mutations.append(
                    {
                        "gene_id": gene_id,
                        "pos": idx + 1,
                        "ref": ref_nt,
                        "alt": alt_nt,
                        "aa_change": aa_change,
                        "mutation_type": mutation_type,
                        "region_type": region_type,
                        "region_name": region_name,
                        "source_sequence": sequence,
                        "target_drug_id": target_drug_id,
                        "status": "PENDING",
                    }
                )
        return mutations

    def bulk_insert_mutations(
        self,
        mutations: Iterable[Dict[str, Any]],
        chunk_size: int = 5000,
        evo_model: str = "evo2_7b",
    ) -> int:
        if self.db is None:
            raise ValueError("DatabaseManager is required for bulk insert.")
        rows = [
            (
                item["gene_id"],
                item["pos"],
                item["ref"],
                item["alt"],
                item.get("aa_change"),
                item.get("region_type", "CDS"),
                item.get("region_name"),
                item.get("source_sequence"),
                item.get("evo_model", evo_model),
                item.get("target_drug_id"),
                item.get("status", "PENDING"),
            )
            for item in mutations
        ]
        if not rows:
            return 0
        before = self.db.conn.total_changes
        sql = "\n            INSERT OR IGNORE INTO mutations\n            (gene_id, pos, ref, alt, aa_change, region_type, region_name, source_sequence, evo_model, target_drug_id, status)\n            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)\n        "
        cursor = self.db.conn.cursor()
        for start in range(0, len(rows), chunk_size):
            chunk = rows[start : start + chunk_size]
            cursor.executemany(sql, chunk)
        self.db.conn.commit()
        return self.db.conn.total_changes - before
