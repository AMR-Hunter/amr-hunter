from __future__ import annotations
import errno
import logging
import os
import shutil
import sqlite3
import threading
import time
from pathlib import Path
from typing import Optional

logger = logging.getLogger("amr_hunter.database")


def _is_disk_quota_error(exc: BaseException) -> bool:
    current: Optional[BaseException] = exc
    quota_errnos = {122, getattr(errno, "EDQUOT", 122)}
    while current is not None:
        if (
            isinstance(current, OSError)
            and getattr(current, "errno", None) in quota_errnos
        ):
            return True
        if "Disk quota exceeded" in str(current):
            return True
        current = current.__cause__ or current.__context__
    return False


class DatabaseManager:
    CORE_JUDGMENT_FIELDS = ("evo_delta", "boltz_ptm", "boltz_iptm", "boltz_plddt")
    SPECIALIZED_SUPPORT_FIELDS = ("boltz_stability_score", "boltz_clash")
    AUXILIARY_ARCHIVE_FIELDS = (
        "boltz_binding_affinity",
        "boltz_affinity_pred_value",
        "boltz_pair_energy",
        "boltz_complex_energy",
        "boltz_artifact_dir",
    )
    MUTATION_STATUSES = (
        "PENDING",
        "EVO_DONE",
        "BOLTZ_READY",
        "LETHAL_SKIP",
        "LOSS_OF_FUNCTION",
        "COMPLETED",
        "RE_SENSITIZED",
        "FAILED",
    )

    def __init__(self, local_db_path: str, share_db_path: str) -> None:
        self.local_db_path = Path(local_db_path)
        self.share_db_path = Path(share_db_path)
        self._share_sync_same_path = self.local_db_path.resolve(
            strict=False
        ) == self.share_db_path.resolve(strict=False)
        self._share_sync_lock = threading.Lock()
        self._share_sync_thread: Optional[threading.Thread] = None
        self._pending_share_snapshot: Optional[Path] = None
        self._share_sync_error: Optional[Exception] = None
        self._share_sync_suspended_reason: Optional[str] = None
        self._bootstrap_local_db()
        self.conn = sqlite3.connect(self.local_db_path)
        self.conn.row_factory = sqlite3.Row
        self._enable_performance_pragmas()
        self._create_schema()
        if self._share_sync_same_path:
            logger.info(
                "[数据库] local_db_path 与 share_db_path 相同；共享同步退化为本地 commit 以避免同路径快照替换覆盖活动连接"
            )

    def _bootstrap_local_db(self) -> None:
        self.local_db_path.parent.mkdir(parents=True, exist_ok=True)
        self.share_db_path.parent.mkdir(parents=True, exist_ok=True)
        if not self.local_db_path.exists() and self.share_db_path.exists():
            shutil.copy2(self.share_db_path, self.local_db_path)

    def _enable_performance_pragmas(self) -> None:
        self.conn.execute("PRAGMA journal_mode=WAL;")
        self.conn.execute("PRAGMA synchronous=NORMAL;")
        self.conn.execute("PRAGMA foreign_keys=ON;")

    def _create_schema(self) -> None:
        self.conn.executescript(
            "\n            CREATE TABLE IF NOT EXISTS species (\n                id INTEGER PRIMARY KEY AUTOINCREMENT,\n                name TEXT NOT NULL UNIQUE,\n                genome_path TEXT\n            );\n\n            CREATE TABLE IF NOT EXISTS drugs (\n                id INTEGER PRIMARY KEY AUTOINCREMENT,\n                name TEXT NOT NULL UNIQUE,\n                sdf_path TEXT\n            );\n\n            CREATE TABLE IF NOT EXISTS genes (\n                id INTEGER PRIMARY KEY AUTOINCREMENT,\n                species_id INTEGER NOT NULL,\n                name TEXT NOT NULL,\n                sequence TEXT,\n                wt_evo_score REAL,\n                pdb_path TEXT,\n                logic_type TEXT DEFAULT 'POSITIVE' CHECK (logic_type IN ('POSITIVE', 'NEGATIVE', 'STRUCTURAL')),\n                scenario TEXT CHECK (scenario IN ('PROTEIN_LIGAND', 'PROTEIN_PROTEIN', 'PROTEIN_DNA', 'FOLDING_ONLY', 'DIMER_DNA', NULL)),\n                pathway_group TEXT,\n                dna_sequence TEXT,\n                ligand_sdf_path TEXT,\n                FOREIGN KEY(species_id) REFERENCES species(id) ON DELETE CASCADE,\n                UNIQUE(species_id, name)\n            );\n\n            -- 字段分组约定：\n            -- 1) 核心判读字段：evo_delta, boltz_ptm, boltz_iptm, boltz_plddt\n            -- 2) 专用支持字段：boltz_stability_score, boltz_clash\n            -- 3) 辅助留档字段：boltz_binding_affinity, boltz_affinity_pred_value,\n            --    boltz_pair_energy, boltz_complex_energy, boltz_artifact_dir\n            CREATE TABLE IF NOT EXISTS mutations (\n                id INTEGER PRIMARY KEY AUTOINCREMENT,\n                gene_id INTEGER NOT NULL,\n                pos INTEGER NOT NULL,\n                ref TEXT NOT NULL,\n                alt TEXT NOT NULL,\n                aa_change TEXT,\n                region_type TEXT NOT NULL DEFAULT 'CDS'\n                    CHECK (region_type IN ('CDS', 'REGULATORY')),\n                region_name TEXT,\n                source_sequence TEXT,\n\n                evo_model TEXT NOT NULL DEFAULT 'evo2_7b',\n                evo_mt_score REAL,\n                evo_delta REAL,\n\n                boltz_binding_affinity REAL,\n                boltz_stability_score REAL,\n                boltz_clash REAL,\n                boltz_iptm REAL,\n                boltz_plddt REAL,\n                boltz_complex_energy REAL,\n                boltz_ptm REAL,\n                boltz_affinity_pred_value REAL,\n                boltz_pair_energy REAL,\n                boltz_artifact_dir TEXT,\n                synergy_score REAL NOT NULL DEFAULT 0,\n                synergy_tags TEXT,\n\n                final_interpretation TEXT,\n                resistance_phenotype TEXT\n                    CHECK (resistance_phenotype IN ('RESISTANT', 'SENSITIVE', 'UNKNOWN') OR resistance_phenotype IS NULL),\n                functional_state TEXT\n                    CHECK (functional_state IN ('INTACT', 'LOSS_OF_FUNCTION', 'LETHAL', 'RE_SENSITIZED', 'UNKNOWN') OR functional_state IS NULL),\n                evidence_code TEXT,\n                target_drug_id INTEGER,\n                retry_count INTEGER NOT NULL DEFAULT 0,\n                last_error TEXT,\n\n                status TEXT NOT NULL DEFAULT 'PENDING'\n                    CHECK (status IN ('PENDING', 'EVO_DONE', 'BOLTZ_READY', 'LETHAL_SKIP', 'LOSS_OF_FUNCTION', 'COMPLETED', 'RE_SENSITIZED', 'FAILED')),\n\n                created_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,\n                updated_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,\n\n                FOREIGN KEY(gene_id) REFERENCES genes(id) ON DELETE CASCADE,\n                FOREIGN KEY(target_drug_id) REFERENCES drugs(id) ON DELETE SET NULL,\n                UNIQUE(gene_id, region_type, region_name, pos, ref, alt, target_drug_id, evo_model)\n            );\n\n            -- Boltz2 结构预测结果缓存，按 (蛋白序列哈希, 场景, 配体哈希) 去重，\n            -- 防止相同结构预测因不同 evo 模型被重复执行。\n            CREATE TABLE IF NOT EXISTS boltz_cache (\n                id INTEGER PRIMARY KEY AUTOINCREMENT,\n                seq_hash TEXT NOT NULL,\n                scenario TEXT NOT NULL,\n                ligand_hash TEXT NOT NULL DEFAULT '',\n                iptm REAL,\n                plddt REAL,\n                complex_energy REAL,\n                binding_affinity REAL,\n                stability_score REAL,\n                clash_score REAL,\n                cached_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,\n                UNIQUE(seq_hash, scenario, ligand_hash)\n            );\n\n            -- 协同/上位效应事件表：记录每次规则命中（多基因协同 或 通路校正）\n            CREATE TABLE IF NOT EXISTS synergy_events (\n                id INTEGER PRIMARY KEY AUTOINCREMENT,\n                rule_name TEXT NOT NULL,\n                rule_type TEXT NOT NULL,\n                rule_source TEXT NOT NULL DEFAULT 'synergy_rules'\n                    CHECK (rule_source IN ('synergy_rules', 'pathway_rules')),\n                description TEXT,\n                synergy_score_modifier REAL NOT NULL DEFAULT 0,\n                matched_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP\n            );\n\n            -- 事件-突变关联表：精确记录哪些突变参与了哪个协同/上位事件\n            CREATE TABLE IF NOT EXISTS synergy_event_mutations (\n                id INTEGER PRIMARY KEY AUTOINCREMENT,\n                event_id INTEGER NOT NULL,\n                mutation_id INTEGER NOT NULL,\n                role TEXT NOT NULL CHECK (role IN ('causal', 'affected')),\n                status_at_match TEXT,\n                FOREIGN KEY(event_id) REFERENCES synergy_events(id) ON DELETE CASCADE,\n                FOREIGN KEY(mutation_id) REFERENCES mutations(id) ON DELETE CASCADE,\n                UNIQUE(event_id, mutation_id)\n            );\n            "
        )
        self._migrate_schema_if_needed()
        self.conn.commit()

    def _migrate_schema_if_needed(self) -> None:
        gene_columns = {
            row["name"]
            for row in self.conn.execute("PRAGMA table_info(genes)").fetchall()
        }
        if "scenario" not in gene_columns:
            self.conn.execute(
                "ALTER TABLE genes ADD COLUMN scenario TEXT CHECK (scenario IN ('PROTEIN_LIGAND', 'PROTEIN_PROTEIN', 'PROTEIN_DNA', 'FOLDING_ONLY', 'DIMER_DNA', NULL))"
            )
        if "pathway_group" not in gene_columns:
            self.conn.execute("ALTER TABLE genes ADD COLUMN pathway_group TEXT")
        if "dna_sequence" not in gene_columns:
            self.conn.execute("ALTER TABLE genes ADD COLUMN dna_sequence TEXT")
        if "ligand_sdf_path" not in gene_columns:
            self.conn.execute("ALTER TABLE genes ADD COLUMN ligand_sdf_path TEXT")
        self.conn.execute(
            "UPDATE genes SET logic_type = 'POSITIVE' WHERE logic_type = 'FORWARD'"
        )
        self.conn.execute(
            "UPDATE genes SET logic_type = 'NEGATIVE' WHERE logic_type = 'REVERSE'"
        )
        genes_table_sql = self.conn.execute(
            "SELECT sql FROM sqlite_master WHERE type='table' AND name='genes'"
        ).fetchone()
        genes_sql_text = genes_table_sql["sql"] or "" if genes_table_sql else ""
        needs_genes_rebuild = (
            "STRUCTURAL" not in genes_sql_text
            or "FOLDING_ONLY" not in genes_sql_text
            or "DIMER_DNA" not in genes_sql_text
        )
        if needs_genes_rebuild:
            self.conn.execute("PRAGMA foreign_keys=OFF")
            self.conn.executescript(
                "\n                CREATE TABLE genes_v21 (\n                    id INTEGER PRIMARY KEY AUTOINCREMENT,\n                    species_id INTEGER NOT NULL,\n                    name TEXT NOT NULL,\n                    sequence TEXT,\n                    wt_evo_score REAL,\n                    pdb_path TEXT,\n                    logic_type TEXT DEFAULT 'POSITIVE' CHECK (logic_type IN ('POSITIVE', 'NEGATIVE', 'STRUCTURAL')),\n                    scenario TEXT CHECK (scenario IN ('PROTEIN_LIGAND', 'PROTEIN_PROTEIN', 'PROTEIN_DNA', 'FOLDING_ONLY', 'DIMER_DNA', NULL)),\n                    pathway_group TEXT,\n                    dna_sequence TEXT,\n                    ligand_sdf_path TEXT,\n                    FOREIGN KEY(species_id) REFERENCES species(id) ON DELETE CASCADE,\n                    UNIQUE(species_id, name)\n                );\n\n                INSERT INTO genes_v21 (\n                    id, species_id, name, sequence, wt_evo_score, pdb_path,\n                    logic_type, scenario, pathway_group, dna_sequence, ligand_sdf_path\n                )\n                SELECT\n                    id, species_id, name, sequence, wt_evo_score, pdb_path,\n                    logic_type, scenario, pathway_group, dna_sequence, ligand_sdf_path\n                FROM genes;\n\n                DROP TABLE genes;\n                ALTER TABLE genes_v21 RENAME TO genes;\n                "
            )
            self.conn.execute("PRAGMA foreign_keys=ON")
        columns = {
            row["name"]
            for row in self.conn.execute("PRAGMA table_info(mutations)").fetchall()
        }
        if "final_interpretation" not in columns:
            self.conn.execute(
                "ALTER TABLE mutations ADD COLUMN final_interpretation TEXT"
            )
        if "boltz_iptm" not in columns:
            self.conn.execute("ALTER TABLE mutations ADD COLUMN boltz_iptm REAL")
        if "boltz_plddt" not in columns:
            self.conn.execute("ALTER TABLE mutations ADD COLUMN boltz_plddt REAL")
        if "boltz_complex_energy" not in columns:
            self.conn.execute(
                "ALTER TABLE mutations ADD COLUMN boltz_complex_energy REAL"
            )
        if "boltz_ptm" not in columns:
            self.conn.execute("ALTER TABLE mutations ADD COLUMN boltz_ptm REAL")
        if "boltz_affinity_pred_value" not in columns:
            self.conn.execute(
                "ALTER TABLE mutations ADD COLUMN boltz_affinity_pred_value REAL"
            )
        if "boltz_pair_energy" not in columns:
            self.conn.execute("ALTER TABLE mutations ADD COLUMN boltz_pair_energy REAL")
        if "boltz_artifact_dir" not in columns:
            self.conn.execute(
                "ALTER TABLE mutations ADD COLUMN boltz_artifact_dir TEXT"
            )
        if "retry_count" not in columns:
            self.conn.execute(
                "ALTER TABLE mutations ADD COLUMN retry_count INTEGER NOT NULL DEFAULT 0"
            )
        if "last_error" not in columns:
            self.conn.execute("ALTER TABLE mutations ADD COLUMN last_error TEXT")
        if "region_type" not in columns:
            self.conn.execute(
                "ALTER TABLE mutations ADD COLUMN region_type TEXT NOT NULL DEFAULT 'CDS' CHECK (region_type IN ('CDS', 'REGULATORY'))"
            )
        if "region_name" not in columns:
            self.conn.execute("ALTER TABLE mutations ADD COLUMN region_name TEXT")
        if "source_sequence" not in columns:
            self.conn.execute("ALTER TABLE mutations ADD COLUMN source_sequence TEXT")
        if "synergy_score" not in columns:
            self.conn.execute(
                "ALTER TABLE mutations ADD COLUMN synergy_score REAL NOT NULL DEFAULT 0"
            )
        if "synergy_tags" not in columns:
            self.conn.execute("ALTER TABLE mutations ADD COLUMN synergy_tags TEXT")
        if "resistance_phenotype" not in columns:
            self.conn.execute(
                "ALTER TABLE mutations ADD COLUMN resistance_phenotype TEXT CHECK (resistance_phenotype IN ('RESISTANT', 'SENSITIVE', 'UNKNOWN') OR resistance_phenotype IS NULL)"
            )
        if "functional_state" not in columns:
            self.conn.execute(
                "ALTER TABLE mutations ADD COLUMN functional_state TEXT CHECK (functional_state IN ('INTACT', 'LOSS_OF_FUNCTION', 'LETHAL', 'RE_SENSITIZED', 'UNKNOWN') OR functional_state IS NULL)"
            )
        if "evidence_code" not in columns:
            self.conn.execute("ALTER TABLE mutations ADD COLUMN evidence_code TEXT")
        if "evo_model" not in columns:
            self.conn.execute(
                "ALTER TABLE mutations ADD COLUMN evo_model TEXT NOT NULL DEFAULT 'evo2_7b'"
            )
        self.conn.execute(
            "UPDATE mutations SET status = 'LOSS_OF_FUNCTION' WHERE status = 'RESISTANT'"
        )
        self.conn.execute(
            "UPDATE mutations SET status = 'BOLTZ_READY' WHERE status = 'READY_FOR_BOLTZ'"
        )
        self.conn.execute(
            "UPDATE mutations SET status = 'BOLTZ_READY' WHERE status = 'BOLTZ_DEFERRED'"
        )
        mutations_table_sql = self.conn.execute(
            "SELECT sql FROM sqlite_master WHERE type='table' AND name='mutations'"
        ).fetchone()
        mutations_sql_text = (
            mutations_table_sql["sql"] or "" if mutations_table_sql else ""
        )
        needs_mutations_rebuild = (
            "region_type" not in mutations_sql_text
            or "region_name" not in mutations_sql_text
            or "source_sequence" not in mutations_sql_text
            or ("synergy_score" not in mutations_sql_text)
            or ("synergy_tags" not in mutations_sql_text)
            or ("resistance_phenotype" not in mutations_sql_text)
            or ("functional_state" not in mutations_sql_text)
            or ("evidence_code" not in mutations_sql_text)
            or ("evo_model" not in mutations_sql_text)
            or ("boltz_ptm" not in mutations_sql_text)
            or ("boltz_affinity_pred_value" not in mutations_sql_text)
            or ("boltz_pair_energy" not in mutations_sql_text)
            or ("boltz_artifact_dir" not in mutations_sql_text)
            or ("BOLTZ_DEFERRED" in mutations_sql_text)
            or (
                "UNIQUE(gene_id, region_type, region_name, pos, ref, alt, target_drug_id, evo_model)"
                not in mutations_sql_text
            )
        )
        if needs_mutations_rebuild:
            self.conn.execute("PRAGMA foreign_keys=OFF")
            self.conn.executescript(
                "\n                CREATE TABLE mutations_v26 (\n                    id INTEGER PRIMARY KEY AUTOINCREMENT,\n                    gene_id INTEGER NOT NULL,\n                    pos INTEGER NOT NULL,\n                    ref TEXT NOT NULL,\n                    alt TEXT NOT NULL,\n                    aa_change TEXT,\n                    region_type TEXT NOT NULL DEFAULT 'CDS'\n                        CHECK (region_type IN ('CDS', 'REGULATORY')),\n                    region_name TEXT,\n                    source_sequence TEXT,\n\n                    evo_model TEXT NOT NULL DEFAULT 'evo2_7b',\n                    evo_mt_score REAL,\n                    evo_delta REAL,\n\n                    boltz_binding_affinity REAL,\n                    boltz_stability_score REAL,\n                    boltz_clash REAL,\n                    boltz_iptm REAL,\n                    boltz_plddt REAL,\n                    boltz_complex_energy REAL,\n                    boltz_ptm REAL,\n                    boltz_affinity_pred_value REAL,\n                    boltz_pair_energy REAL,\n                    boltz_artifact_dir TEXT,\n                    synergy_score REAL NOT NULL DEFAULT 0,\n                    synergy_tags TEXT,\n\n                    final_interpretation TEXT,\n                    resistance_phenotype TEXT\n                        CHECK (resistance_phenotype IN ('RESISTANT', 'SENSITIVE', 'UNKNOWN') OR resistance_phenotype IS NULL),\n                    functional_state TEXT\n                        CHECK (functional_state IN ('INTACT', 'LOSS_OF_FUNCTION', 'LETHAL', 'RE_SENSITIZED', 'UNKNOWN') OR functional_state IS NULL),\n                    evidence_code TEXT,\n                    target_drug_id INTEGER,\n                    retry_count INTEGER NOT NULL DEFAULT 0,\n                    last_error TEXT,\n\n                    status TEXT NOT NULL DEFAULT 'PENDING'\n                        CHECK (status IN ('PENDING', 'EVO_DONE', 'BOLTZ_READY', 'LETHAL_SKIP', 'LOSS_OF_FUNCTION', 'COMPLETED', 'RE_SENSITIZED', 'FAILED')),\n\n                    created_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,\n                    updated_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,\n\n                    FOREIGN KEY(gene_id) REFERENCES genes(id) ON DELETE CASCADE,\n                    FOREIGN KEY(target_drug_id) REFERENCES drugs(id) ON DELETE SET NULL,\n                    UNIQUE(gene_id, region_type, region_name, pos, ref, alt, target_drug_id, evo_model)\n                );\n\n                INSERT OR IGNORE INTO mutations_v26 (\n                    id, gene_id, pos, ref, alt, aa_change, region_type, region_name, source_sequence,\n                    evo_model, evo_mt_score, evo_delta,\n                    boltz_binding_affinity, boltz_stability_score, boltz_clash,\n                    boltz_iptm, boltz_plddt, boltz_complex_energy,\n                    boltz_ptm, boltz_affinity_pred_value, boltz_pair_energy,\n                    boltz_artifact_dir,\n                    synergy_score, synergy_tags,\n                    final_interpretation, resistance_phenotype, functional_state, evidence_code,\n                    target_drug_id, retry_count, last_error,\n                    status, created_at, updated_at\n                )\n                SELECT\n                    id,\n                    gene_id,\n                    pos,\n                    ref,\n                    alt,\n                    aa_change,\n                    COALESCE(region_type, 'CDS'),\n                    region_name,\n                    source_sequence,\n                    COALESCE(evo_model, 'evo2_7b'),\n                    evo_mt_score,\n                    evo_delta,\n                    boltz_binding_affinity,\n                    boltz_stability_score,\n                    boltz_clash,\n                    boltz_iptm,\n                    boltz_plddt,\n                    boltz_complex_energy,\n                    boltz_ptm,\n                    boltz_affinity_pred_value,\n                    boltz_pair_energy,\n                    boltz_artifact_dir,\n                    COALESCE(synergy_score, 0),\n                    synergy_tags,\n                    final_interpretation,\n                    resistance_phenotype,\n                    functional_state,\n                    evidence_code,\n                    target_drug_id,\n                    COALESCE(retry_count, 0),\n                    last_error,\n                    status,\n                    created_at,\n                    updated_at\n                FROM mutations;\n\n                DROP TABLE mutations;\n                ALTER TABLE mutations_v26 RENAME TO mutations;\n                "
            )
            self.conn.execute("PRAGMA foreign_keys=ON")
        self.conn.execute(
            "\n            CREATE TABLE IF NOT EXISTS boltz_cache (\n                id INTEGER PRIMARY KEY AUTOINCREMENT,\n                seq_hash TEXT NOT NULL,\n                scenario TEXT NOT NULL,\n                ligand_hash TEXT NOT NULL DEFAULT '',\n                iptm REAL,\n                plddt REAL,\n                complex_energy REAL,\n                binding_affinity REAL,\n                stability_score REAL,\n                ptm REAL,\n                affinity_pred_value REAL,\n                pair_energy REAL,\n                clash_score REAL,\n                cached_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,\n                UNIQUE(seq_hash, scenario, ligand_hash)\n            )\n            "
        )
        boltz_cache_columns = {
            row["name"]
            for row in self.conn.execute("PRAGMA table_info(boltz_cache)").fetchall()
        }
        if "ptm" not in boltz_cache_columns:
            self.conn.execute("ALTER TABLE boltz_cache ADD COLUMN ptm REAL")
        if "affinity_pred_value" not in boltz_cache_columns:
            self.conn.execute(
                "ALTER TABLE boltz_cache ADD COLUMN affinity_pred_value REAL"
            )
        if "pair_energy" not in boltz_cache_columns:
            self.conn.execute("ALTER TABLE boltz_cache ADD COLUMN pair_energy REAL")
        self.conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_genes_species_id ON genes(species_id)"
        )
        self.conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_genes_pathway_group ON genes(pathway_group)"
        )
        self.conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_mutations_gene_id ON mutations(gene_id)"
        )
        self.conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_mutations_target_drug_id ON mutations(target_drug_id)"
        )
        self.conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_mutations_status ON mutations(status)"
        )
        self.conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_mutations_evo_model ON mutations(evo_model)"
        )
        self.conn.commit()

    def clear_synergy_events(self) -> None:
        self.conn.execute("DELETE FROM synergy_event_mutations")
        self.conn.execute("DELETE FROM synergy_events")
        self.conn.commit()

    def _raise_share_sync_error_if_any(self) -> None:
        if self._share_sync_error is not None:
            raise RuntimeError(
                f"共享盘数据库异步同步失败: {self._share_sync_error}"
            ) from self._share_sync_error

    def _store_unsynced_local_snapshot(self, snapshot_path: Path) -> Optional[Path]:
        fallback_path = self.local_db_path.with_name(
            f"{self.local_db_path.stem}.share_unsynced{self.local_db_path.suffix}"
        )
        try:
            shutil.copy2(snapshot_path, fallback_path)
            return fallback_path
        except Exception as exc:
            logger.warning("[数据库] 保存本地未同步快照失败: %s", exc)
            return None

    def _create_share_snapshot(self) -> Path:
        self.conn.commit()
        self.conn.execute("PRAGMA wal_checkpoint(TRUNCATE)")
        self.conn.commit()
        snapshot_path = self.local_db_path.with_name(
            f"{self.local_db_path.stem}.share_sync_snapshot.{time.monotonic_ns()}{self.local_db_path.suffix}"
        )
        with sqlite3.connect(self.local_db_path) as src_conn:
            with sqlite3.connect(snapshot_path) as snapshot_conn:
                src_conn.backup(snapshot_conn)
        return snapshot_path

    def _copy_snapshot_to_share(self, snapshot_path: Path) -> None:
        self.share_db_path.parent.mkdir(parents=True, exist_ok=True)
        tmp_share_path = self.share_db_path.with_suffix(".db.tmp")
        tmp_share_path.unlink(missing_ok=True)
        try:
            shutil.copy2(snapshot_path, tmp_share_path)
            os.replace(tmp_share_path, self.share_db_path)
        finally:
            tmp_share_path.unlink(missing_ok=True)

    def _share_sync_worker(self, snapshot_path: Path) -> None:
        current_snapshot: Optional[Path] = snapshot_path
        while current_snapshot is not None:
            next_snapshot: Optional[Path] = None
            try:
                self._copy_snapshot_to_share(current_snapshot)
            except Exception as exc:
                if _is_disk_quota_error(exc):
                    fallback_path = self._store_unsynced_local_snapshot(
                        current_snapshot
                    )
                    with self._share_sync_lock:
                        self._share_sync_suspended_reason = str(exc)
                        next_snapshot = self._pending_share_snapshot
                        self._pending_share_snapshot = None
                        self._share_sync_thread = None
                    if next_snapshot is not None:
                        next_snapshot.unlink(missing_ok=True)
                    logger.warning(
                        "[数据库] 共享盘配额不足，数据库共享同步已暂停；本地未同步快照保留在 %s",
                        fallback_path or current_snapshot,
                    )
                    return
                with self._share_sync_lock:
                    self._share_sync_error = exc
                    next_snapshot = self._pending_share_snapshot
                    self._pending_share_snapshot = None
                    self._share_sync_thread = None
                if next_snapshot is not None:
                    next_snapshot.unlink(missing_ok=True)
                return
            finally:
                current_snapshot.unlink(missing_ok=True)
            with self._share_sync_lock:
                next_snapshot = self._pending_share_snapshot
                self._pending_share_snapshot = None
                if next_snapshot is None:
                    self._share_sync_thread = None
            current_snapshot = next_snapshot

    def sync_to_share_async(self) -> None:
        if self._share_sync_same_path:
            self.conn.commit()
            return
        if self._share_sync_suspended_reason is not None:
            return
        self._raise_share_sync_error_if_any()
        snapshot_path = self._create_share_snapshot()
        with self._share_sync_lock:
            self._raise_share_sync_error_if_any()
            if (
                self._share_sync_thread is not None
                and self._share_sync_thread.is_alive()
            ):
                previous_snapshot = self._pending_share_snapshot
                self._pending_share_snapshot = snapshot_path
                if previous_snapshot is not None:
                    previous_snapshot.unlink(missing_ok=True)
                return
            self._share_sync_thread = threading.Thread(
                target=self._share_sync_worker,
                args=(snapshot_path,),
                name="db-share-sync",
                daemon=True,
            )
            self._share_sync_thread.start()

    def wait_for_pending_share_sync(self) -> None:
        while True:
            with self._share_sync_lock:
                thread = self._share_sync_thread
            if thread is None:
                break
            thread.join()
        self._raise_share_sync_error_if_any()

    def sync_to_share(self) -> None:
        self.sync_to_share_async()
        self.wait_for_pending_share_sync()

    def close(self) -> None:
        if self.conn:
            self.wait_for_pending_share_sync()
            self.conn.commit()
            self.conn.close()

    def __enter__(self) -> "DatabaseManager":
        return self

    def __exit__(self, exc_type, exc, exc_tb) -> None:
        self.close()


def build_default_manager(
    local_db_path: str = "/root/amr_hunter/data/db/amr_hunter.db",
    share_db_path: str = "/root/gpufree-share/amr_hunter/db/amr_hunter.db",
) -> DatabaseManager:
    return DatabaseManager(local_db_path=local_db_path, share_db_path=share_db_path)
