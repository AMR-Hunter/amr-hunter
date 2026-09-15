from __future__ import annotations
import argparse
import asyncio
import contextlib
from dataclasses import asdict, dataclass
import errno
import hashlib
import io
import json
import logging
import os
from queue import Queue
import shutil
import subprocess
import sys
import shlex
import tarfile
import threading
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

try:
    from fastapi import Request
except ImportError:
    Request = object
logger = logging.getLogger("boltz2_server")
logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
_model_dir: Path = Path("/root/gpufree-data/models/boltz2")
_ready: bool = False
_startup_error: Optional[str] = None
_job_lock_guard = threading.Lock()
_job_locks: Dict[str, threading.Lock] = {}
_predict_concurrency = max(1, int(os.getenv("BOLTZ_MAX_CONCURRENT_PREDICTS", "1")))
_predict_slots = threading.BoundedSemaphore(_predict_concurrency)


def _env_flag(name: str, default: bool) -> bool:
    raw_value = os.getenv(name)
    if raw_value is None:
        return default
    return raw_value.strip().lower() not in {"0", "false", "no", "off"}


_prepare_concurrency = max(1, int(os.getenv("BOLTZ_MAX_CONCURRENT_PREPARES", "1")))
_prepare_long_job_threshold = float(
    os.getenv("BOLTZ_PREPARE_LONG_JOB_THRESHOLD", "35.0")
)
_prepare_age_priority_seconds = float(
    os.getenv("BOLTZ_PREPARE_AGE_PRIORITY_SECONDS", "120.0")
)
_prepare_long_job_cap = max(
    1,
    min(
        _prepare_concurrency,
        int(
            os.getenv(
                "BOLTZ_MAX_CONCURRENT_LONG_PREPARES",
                str(max(1, _prepare_concurrency - 1)),
            )
        ),
    ),
)
_prepare_scheduler_guard = threading.Condition()
_prepare_scheduler_sequence = 0
_prepare_active_jobs = 0
_prepare_active_long_jobs = 0
_prepare_waiters: List["_PrepareTicket"] = []
_remote_msa_concurrency = max(1, int(os.getenv("BOLTZ_MAX_CONCURRENT_REMOTE_MSA", "1")))
_remote_msa_slots = threading.BoundedSemaphore(_remote_msa_concurrency)
_remote_msa_guard = threading.Lock()
_remote_msa_active_jobs = 0
_persistent_runtime_enabled = _env_flag("BOLTZ_PERSISTENT_RUNTIME", True)
_runtime_pool_guard = threading.Lock()
_runtime_pool: Optional[Queue["_PersistentBoltzRunner"]] = None
_runtime_pool_ready = False
_runtime_init_error: Optional[str] = None
_runtime_runner_count = 0


@dataclass
class _PrepareTicket:
    sequence: int
    enqueued_at: float
    job_key: str
    duration_score: float
    is_long: bool


def _payload_duration_score(payload: Dict[str, Any]) -> float:
    meta = payload.get("_meta") if isinstance(payload.get("_meta"), dict) else {}
    raw_value = meta.get("estimated_duration_score") if isinstance(meta, dict) else 0.0
    try:
        return float(raw_value or 0.0)
    except (TypeError, ValueError):
        return 0.0


def _select_prepare_ticket(now: float) -> Optional[_PrepareTicket]:
    if _prepare_active_jobs >= _prepare_concurrency or not _prepare_waiters:
        return None
    short_waiters = [ticket for ticket in _prepare_waiters if not ticket.is_long]
    eligible_waiters = list(_prepare_waiters)
    if _prepare_active_long_jobs >= _prepare_long_job_cap:
        if short_waiters:
            eligible_waiters = short_waiters
    aged_waiters = [
        ticket
        for ticket in eligible_waiters
        if now - ticket.enqueued_at >= _prepare_age_priority_seconds
    ]
    if aged_waiters:
        return min(
            aged_waiters, key=lambda ticket: (ticket.enqueued_at, ticket.sequence)
        )
    return min(
        eligible_waiters,
        key=lambda ticket: (
            1 if ticket.is_long else 0,
            ticket.duration_score,
            ticket.enqueued_at,
            ticket.sequence,
        ),
    )


@contextlib.contextmanager
def _acquire_prepare_slot(job_key: str, duration_score: float):
    global _prepare_scheduler_sequence, _prepare_active_jobs, _prepare_active_long_jobs
    with _prepare_scheduler_guard:
        ticket = _PrepareTicket(
            sequence=_prepare_scheduler_sequence,
            enqueued_at=time.perf_counter(),
            job_key=job_key,
            duration_score=duration_score,
            is_long=duration_score >= _prepare_long_job_threshold,
        )
        _prepare_scheduler_sequence += 1
        _prepare_waiters.append(ticket)
        _prepare_scheduler_guard.notify_all()
        while True:
            now = time.perf_counter()
            selected_ticket = _select_prepare_ticket(now)
            if selected_ticket is ticket:
                _prepare_waiters.remove(ticket)
                _prepare_active_jobs += 1
                if ticket.is_long:
                    _prepare_active_long_jobs += 1
                queue_wait = now - ticket.enqueued_at
                logger.info(
                    "[Boltz2] 进入输入处理槽位: job=%s queue_wait=%.2fs dur_score=%.2f long=%s active=%d/%d active_long=%d/%d pending=%d",
                    ticket.job_key[:12],
                    queue_wait,
                    ticket.duration_score,
                    ticket.is_long,
                    _prepare_active_jobs,
                    _prepare_concurrency,
                    _prepare_active_long_jobs,
                    _prepare_long_job_cap,
                    len(_prepare_waiters),
                )
                break
            _prepare_scheduler_guard.wait(timeout=1.0)
    try:
        yield queue_wait
    finally:
        with _prepare_scheduler_guard:
            _prepare_active_jobs = max(0, _prepare_active_jobs - 1)
            if ticket.is_long:
                _prepare_active_long_jobs = max(0, _prepare_active_long_jobs - 1)
            logger.info(
                "[Boltz2] 释放输入处理槽位: job=%s active=%d/%d active_long=%d/%d pending=%d",
                ticket.job_key[:12],
                _prepare_active_jobs,
                _prepare_concurrency,
                _prepare_active_long_jobs,
                _prepare_long_job_cap,
                len(_prepare_waiters),
            )
            _prepare_scheduler_guard.notify_all()


class BoltzPredictError(RuntimeError):

    def __init__(self, message: str, *, stdout: str = "", stderr: str = "") -> None:
        super().__init__(message)
        self.stdout = stdout
        self.stderr = stderr


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


def _check_boltz_installed() -> bool:
    try:
        import boltz

        return True
    except ImportError:
        return False


def _prefetch_weights(model_dir: Path) -> None:
    global _ready, _startup_error
    if not _check_boltz_installed():
        _startup_error = "boltz 包未安装，请运行：pip install boltz"
        logger.error(_startup_error)
        return
    boltz_cache = Path(os.getenv("BOLTZ_CACHE_DIR", str(model_dir)))
    boltz_home = Path.home() / ".boltz"
    if model_dir.exists() and any(model_dir.iterdir()):
        if not boltz_home.exists() or not any(boltz_home.iterdir()):
            logger.info(
                "[Boltz2] 拷贝权重到 boltz 默认缓存目录 %s → %s", model_dir, boltz_home
            )
            boltz_home.parent.mkdir(parents=True, exist_ok=True)
            if boltz_home.exists():
                shutil.rmtree(boltz_home)
            shutil.copytree(model_dir, boltz_home)
            logger.info("[Boltz2] 权重拷贝完成")
        else:
            logger.info("[Boltz2] ~/.boltz 已有权重，跳过拷贝")
    else:
        logger.info("[Boltz2] model_dir 为空，将由 boltz 在首次运行时自动下载权重")
    try:
        import boltz.main

        logger.info("[Boltz2] boltz 包验证通过")
        _ready = True
    except Exception as exc:
        _startup_error = f"boltz 包验证失败：{exc}"
        logger.error(_startup_error)


def _normalize_protein_sequence(sequence: Any) -> str:
    return "".join(str(sequence or "").split()).upper()


def _resolve_explicit_msa_path(protein_info: Dict[str, Any]) -> Optional[Path]:
    for key in ("msa", "msa_path"):
        raw_value = protein_info.get(key)
        if not raw_value:
            continue
        candidate = Path(str(raw_value)).expanduser()
        if candidate.exists():
            return candidate
        logger.warning("[Boltz2] 指定的 MSA 文件不存在：%s", candidate)
    return None


def _sequence_msa_cache_enabled() -> bool:
    return _env_flag("BOLTZ_SEQUENCE_MSA_CACHE", True)


def _local_msa_enabled() -> bool:
    return _env_flag("BOLTZ_LOCAL_MSA_ENABLED", False)


def _local_msa_require_local() -> bool:
    return _env_flag("BOLTZ_LOCAL_MSA_REQUIRE_LOCAL", False)


def _get_local_msa_search_bin() -> str:
    return (
        os.getenv("BOLTZ_LOCAL_MSA_BIN", "colabfold_search").strip()
        or "colabfold_search"
    )


def _get_local_msa_db_dir() -> Optional[Path]:
    configured = os.getenv("BOLTZ_LOCAL_MSA_DB_DIR", "").strip()
    if not configured:
        return None
    return Path(configured).expanduser()


def _get_local_msa_mmseqs_bin() -> Optional[str]:
    configured = os.getenv("BOLTZ_LOCAL_MSA_MMSEQS_BIN", "").strip()
    return configured or None


def _get_local_msa_extra_args() -> List[str]:
    configured = os.getenv("BOLTZ_LOCAL_MSA_EXTRA_ARGS", "").strip()
    return shlex.split(configured) if configured else []


def _get_local_msa_keep_tmp() -> bool:
    return _env_flag("BOLTZ_LOCAL_MSA_KEEP_TMP", False)


def _resolve_local_msa_work_dir(base_dir: Path, cache_key: str) -> Path:
    configured = os.getenv("BOLTZ_LOCAL_MSA_WORK_DIR", "").strip()
    root = Path(configured).expanduser() if configured else base_dir / "__local_msa"
    work_dir = root / cache_key
    work_dir.mkdir(parents=True, exist_ok=True)
    return work_dir


def _local_msa_failure(message: str) -> Optional[Path]:
    if _local_msa_require_local():
        raise RuntimeError(message)
    logger.warning("[Boltz2] 跳过本地 MSA 生成：%s", message)
    return None


def _command_available(command_text: str) -> bool:
    candidate = Path(command_text).expanduser()
    return candidate.exists() or shutil.which(command_text) is not None


def _get_boltz_stable_msa_cache_root() -> Path:
    configured = os.getenv("BOLTZ_STABLE_MSA_CACHE_DIR", "").strip()
    if configured:
        root = Path(configured).expanduser()
    else:
        root = Path("/root/gpufree-data/amr_hunter_runtime/cache/boltz_sequence_msas")
    root.mkdir(parents=True, exist_ok=True)
    return root


def _protein_sequence_cache_key(sequence: Any) -> str:
    normalized = _normalize_protein_sequence(sequence)
    return hashlib.sha1(normalized.encode("utf-8")).hexdigest()


def _find_cached_sequence_msa(sequence: Any) -> Optional[Path]:
    if not _sequence_msa_cache_enabled():
        return None
    cache_root = _get_boltz_stable_msa_cache_root()
    cache_key = _protein_sequence_cache_key(sequence)
    for suffix in (".csv", ".a3m"):
        candidate = cache_root / f"{cache_key}{suffix}"
        if candidate.exists():
            return candidate
    return None


def _generate_local_sequence_msa(sequence: str, work_root: Path) -> Optional[Path]:
    cached = _find_cached_sequence_msa(sequence)
    if cached is not None:
        return cached
    if not _local_msa_enabled():
        return None
    search_bin = _get_local_msa_search_bin()
    if not _command_available(search_bin):
        return _local_msa_failure(f"本地 MSA 命令不可用：{search_bin}")
    db_dir = _get_local_msa_db_dir()
    if db_dir is None:
        return _local_msa_failure("未设置 BOLTZ_LOCAL_MSA_DB_DIR，本地 MSA 无法运行")
    if not db_dir.exists():
        return _local_msa_failure(f"本地 MSA 数据库目录不存在：{db_dir}")
    cache_key = _protein_sequence_cache_key(sequence)
    work_dir = _resolve_local_msa_work_dir(work_root, cache_key)
    query_path = work_dir / f"{cache_key}.fasta"
    output_dir = work_dir / "out"
    if output_dir.exists():
        shutil.rmtree(output_dir, ignore_errors=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    query_path.write_text(f">{cache_key}\n{sequence}\n", encoding="utf-8")
    command = [search_bin]
    mmseqs_bin = _get_local_msa_mmseqs_bin()
    if mmseqs_bin:
        command.extend(["--mmseqs", mmseqs_bin])
    command.extend(_get_local_msa_extra_args())
    command.extend([str(query_path), str(db_dir), str(output_dir)])
    started_at = time.perf_counter()
    try:
        completed = subprocess.run(
            command, capture_output=True, text=True, cwd=str(work_dir), check=False
        )
    except Exception as exc:
        return _local_msa_failure(f"执行本地 MSA 命令失败：{exc}")
    if completed.returncode != 0:
        stdout_tail = completed.stdout[-1200:] if completed.stdout else ""
        stderr_tail = completed.stderr[-1200:] if completed.stderr else ""
        return _local_msa_failure(
            f"本地 MSA 命令返回非零退出码：rc={completed.returncode}\nstdout:\n{stdout_tail}\nstderr:\n{stderr_tail}"
        )
    preferred = output_dir / f"{cache_key}.a3m"
    candidates = []
    if preferred.exists():
        candidates.append(preferred)
    candidates.extend(
        sorted(
            [path for path in output_dir.rglob("*.a3m") if path != preferred],
            key=lambda path: (-path.stat().st_size, path.as_posix()),
        )
    )
    if not candidates:
        return _local_msa_failure(f"本地 MSA 命令未生成任何 a3m 文件：{output_dir}")
    destination = _get_boltz_stable_msa_cache_root() / f"{cache_key}.a3m"
    if not destination.exists():
        shutil.copy2(candidates[0], destination)
        logger.info(
            "[Boltz2] 本地 MSA 生成完成: seq=%s total=%.2fs path=%s",
            cache_key[:12],
            time.perf_counter() - started_at,
            destination,
        )
    if not _get_local_msa_keep_tmp():
        shutil.rmtree(work_dir, ignore_errors=True)
    return destination


def _collect_unique_protein_entities(
    boltz_input: Dict[str, Any],
) -> List[Dict[str, Any]]:
    entities: List[Dict[str, Any]] = []
    seen_sequences = set()
    for item in boltz_input.get("sequences", []):
        protein_info = item.get("protein") if isinstance(item, dict) else None
        if not isinstance(protein_info, dict):
            continue
        sequence = _normalize_protein_sequence(protein_info.get("sequence", ""))
        if not sequence or sequence in seen_sequences:
            continue
        seen_sequences.add(sequence)
        msa_value = protein_info.get("msa")
        msa_path = None
        if msa_value not in (None, "", 0, -1):
            msa_path = Path(str(msa_value)).expanduser()
        entities.append({"sequence": sequence, "msa_path": msa_path})
    return entities


def _resolve_protein_msa_candidate(
    protein_info: Dict[str, Any], *, attach_cached_msas: bool
) -> Optional[Path]:
    explicit_msa_path = _resolve_explicit_msa_path(protein_info)
    if explicit_msa_path is not None:
        return explicit_msa_path
    if not attach_cached_msas:
        return None
    return _find_cached_sequence_msa(protein_info.get("sequence", ""))


def _materialize_local_msas_for_input(
    payload: Dict[str, Any], input_file: Path, work_root: Path
) -> int:
    if not _local_msa_enabled():
        return 0
    base_boltz_input = _build_boltz_input(payload, attach_cached_msas=False)
    generated = 0
    for entity in _collect_unique_protein_entities(base_boltz_input):
        explicit_msa_path = entity.get("msa_path")
        if isinstance(explicit_msa_path, Path) and explicit_msa_path.exists():
            continue
        if _find_cached_sequence_msa(entity["sequence"]) is not None:
            continue
        generated_path = _generate_local_sequence_msa(entity["sequence"], work_root)
        if generated_path is not None and generated_path.exists():
            generated += 1
    if generated:
        input_file.write_text(
            "# boltz input — auto-generated by AMR-Hunter\n"
            + _yaml_dump(_build_boltz_input(payload)),
            encoding="utf-8",
        )
    return generated


def _persist_stable_sequence_msas(
    boltz_input: Dict[str, Any], output_root: Path, input_stem: str
) -> None:
    if not _sequence_msa_cache_enabled():
        return
    cache_root = _get_boltz_stable_msa_cache_root()
    for entity_index, entity in enumerate(
        _collect_unique_protein_entities(boltz_input)
    ):
        cache_key = _protein_sequence_cache_key(entity["sequence"])
        explicit_msa_path = entity.get("msa_path")
        if isinstance(explicit_msa_path, Path) and explicit_msa_path.exists():
            suffix = explicit_msa_path.suffix.lower() or ".a3m"
            destination = cache_root / f"{cache_key}{suffix}"
            if not destination.exists():
                shutil.copy2(explicit_msa_path, destination)
                logger.info(
                    "[Boltz2] 写入稳定 MSA 缓存: seq=%s source=explicit path=%s",
                    cache_key[:12],
                    destination,
                )
            continue
        generated_msa_path = output_root / "msa" / f"{input_stem}_{entity_index}.csv"
        if not generated_msa_path.exists():
            continue
        destination = cache_root / f"{cache_key}.csv"
        if destination.exists():
            continue
        shutil.copy2(generated_msa_path, destination)
        logger.info(
            "[Boltz2] 写入稳定 MSA 缓存: seq=%s source=generated path=%s",
            cache_key[:12],
            destination,
        )


def _boltz_input_needs_remote_msa(boltz_input: Dict[str, Any]) -> bool:
    for entity in _collect_unique_protein_entities(boltz_input):
        msa_path = entity.get("msa_path")
        if not isinstance(msa_path, Path) or not msa_path.exists():
            return True
    return False


@contextlib.contextmanager
def _acquire_remote_msa_slot(job_key: str, *, required: bool):
    global _remote_msa_active_jobs
    if not required:
        yield 0.0
        return
    queued_at = time.perf_counter()
    _remote_msa_slots.acquire()
    queue_wait = time.perf_counter() - queued_at
    with _remote_msa_guard:
        _remote_msa_active_jobs += 1
        active_jobs = _remote_msa_active_jobs
    logger.info(
        "[Boltz2] 进入远端 MSA 槽位: job=%s queue_wait=%.2fs active=%d/%d",
        job_key[:12],
        queue_wait,
        active_jobs,
        _remote_msa_concurrency,
    )
    try:
        yield queue_wait
    finally:
        with _remote_msa_guard:
            _remote_msa_active_jobs = max(0, _remote_msa_active_jobs - 1)
            active_jobs = _remote_msa_active_jobs
        logger.info(
            "[Boltz2] 释放远端 MSA 槽位: job=%s active=%d/%d",
            job_key[:12],
            active_jobs,
            _remote_msa_concurrency,
        )
        _remote_msa_slots.release()


def _build_boltz_input(
    payload: Dict[str, Any], *, attach_cached_msas: bool = True
) -> Dict[str, Any]:
    sequences = []
    properties = []
    scenario = str(payload.get("scenario", "protein_only")).lower()
    if "chains" in payload:
        pending_multichain_msas: List[Optional[Path]] = []
        for i, chain in enumerate(payload["chains"]):
            chain_type = str(chain.get("type", "protein")).lower()
            if chain_type == "protein":
                protein_entry = {"id": chr(65 + i), "sequence": chain["sequence"]}
                pending_multichain_msas.append(
                    _resolve_protein_msa_candidate(
                        chain, attach_cached_msas=attach_cached_msas
                    )
                )
                sequences.append({"protein": protein_entry})
        if pending_multichain_msas:
            if all((path is not None for path in pending_multichain_msas)):
                protein_index = 0
                for item in sequences:
                    protein_info = (
                        item.get("protein") if isinstance(item, dict) else None
                    )
                    if not isinstance(protein_info, dict):
                        continue
                    protein_info["msa"] = str(pending_multichain_msas[protein_index])
                    protein_index += 1
            elif any((path is not None for path in pending_multichain_msas)):
                logger.info(
                    "[Boltz2] 多链输入仅部分蛋白命中显式/缓存 MSA，跳过全部显式 MSA 以避免 custom/auto 混用"
                )
    else:
        protein_info = payload.get("protein", {})
        protein_entry = {"id": "A", "sequence": protein_info.get("sequence", "")}
        msa_path = _resolve_protein_msa_candidate(
            protein_info, attach_cached_msas=attach_cached_msas
        )
        if msa_path is not None:
            protein_entry["msa"] = str(msa_path)
        sequences.append({"protein": protein_entry})
    if scenario == "protein_ligand":
        ligand_info = payload.get("ligand", {})
        sdf_path = ligand_info.get("sdf_path")
        ligand_added = False
        if sdf_path and Path(sdf_path).exists():
            smiles = _sdf_to_smiles(Path(sdf_path))
            if smiles:
                sequences.append({"ligand": {"id": "B", "smiles": smiles}})
                ligand_added = True
            else:
                sequences.append({"ligand": {"id": "B", "ccd": "LIG"}})
                ligand_added = True
        else:
            logger.warning("[Boltz2] 配体 SDF 不存在：%s，跳过配体", sdf_path)
        if ligand_added:
            properties.append({"affinity": {"binder": "B"}})
    elif scenario == "protein_dna":
        dna_info = payload.get("dna", {})
        dna_seq = dna_info.get("sequence", "")
        if dna_seq:
            sequences.append(
                {"dna": {"id": "B", "sequence": dna_seq.upper().replace("U", "T")}}
            )
    boltz_input = {"version": 1, "sequences": sequences}
    if properties:
        boltz_input["properties"] = properties
    template_pdb = None
    if "protein" in payload:
        template_pdb = payload["protein"].get("template_pdb")
    elif "template_pdb" in payload:
        template_pdb = payload["template_pdb"]
    if template_pdb and Path(template_pdb).exists():
        boltz_input["templates"] = [{"pdb": str(template_pdb), "chain_id": "A"}]
    return boltz_input


def _write_boltz_input(payload: Dict[str, Any], work_dir: Path) -> Path:
    boltz_input = _build_boltz_input(payload)
    input_file = work_dir / "input.yaml"
    input_file.write_text(
        "# boltz input — auto-generated by AMR-Hunter\n" + _yaml_dump(boltz_input),
        encoding="utf-8",
    )
    return input_file


def _yaml_dump(data: Any) -> str:
    try:
        import yaml

        return yaml.dump(
            data, allow_unicode=True, sort_keys=False, default_flow_style=False
        )
    except ImportError:
        import json

        return json.dumps(data, ensure_ascii=False, indent=2)


def _sdf_to_smiles(sdf_path: Path) -> Optional[str]:
    try:
        from rdkit import Chem
        from rdkit.Chem.rdmolfiles import MolToSmiles, SDMolSupplier

        suppl = SDMolSupplier(str(sdf_path), removeHs=False)
        mol = next((m for m in suppl if m is not None), None)
        if mol:
            return MolToSmiles(mol)
    except ImportError:
        pass
    except Exception as exc:
        logger.warning("[Boltz2] RDKit SMILES 提取失败：%s", exc)
    try:
        content = sdf_path.read_text(encoding="utf-8", errors="ignore")
        for line in content.splitlines():
            line = line.strip()
            if (
                line
                and (not line.startswith(">"))
                and (not line.startswith("$"))
                and (not line[0].isdigit())
                and (len(line) > 10)
                and any((c in line for c in ("C", "N", "O", "S", "c", "n", "o")))
            ):
                if all((c in "CNOSPFClBrI=#@\\/()+[]0123456789-." for c in line)):
                    return line
    except Exception:
        pass
    return None


def _get_job_lock(job_key: str) -> threading.Lock:
    with _job_lock_guard:
        lock = _job_locks.get(job_key)
        if lock is None:
            lock = threading.Lock()
            _job_locks[job_key] = lock
        return lock


def _gpu_runtime_summary() -> Dict[str, Any]:
    summary: Dict[str, Any] = {
        "torch_version": None,
        "cuda_available": False,
        "device_count": 0,
        "device_name": None,
        "cuda_visible_devices": os.getenv("CUDA_VISIBLE_DEVICES", ""),
    }
    try:
        import torch

        summary["torch_version"] = torch.__version__
        summary["cuda_available"] = bool(torch.cuda.is_available())
        summary["device_count"] = int(torch.cuda.device_count())
        if summary["cuda_available"] and summary["device_count"]:
            summary["device_name"] = torch.cuda.get_device_name(0)
    except Exception as exc:
        summary["error"] = str(exc)
    return summary


def _log_gpu_runtime(context: str) -> Dict[str, Any]:
    summary = _gpu_runtime_summary()
    logger.info(
        "[Boltz2] GPU runtime (%s): cuda_available=%s device_count=%s device=%s visible=%s torch=%s",
        context,
        summary.get("cuda_available"),
        summary.get("device_count"),
        summary.get("device_name") or "<none>",
        summary.get("cuda_visible_devices") or "<unset>",
        summary.get("torch_version") or "<unknown>",
    )
    if summary.get("error"):
        logger.warning(
            "[Boltz2] GPU runtime probe error (%s): %s", context, summary["error"]
        )
    return summary


def _configure_torch_matmul_precision() -> None:
    precision = os.getenv("BOLTZ_FLOAT32_MATMUL_PRECISION", "high").strip().lower()
    if precision not in {"highest", "high", "medium"}:
        logger.warning(
            "[Boltz2] 忽略无效的 float32 matmul precision 配置: %s", precision
        )
        return
    try:
        import torch

        torch.set_float32_matmul_precision(precision)
        logger.info("[Boltz2] torch float32 matmul precision=%s", precision)
    except Exception as exc:
        logger.warning("[Boltz2] 设置 torch float32 matmul precision 失败: %s", exc)


def _resolve_accelerator() -> str:
    allow_cpu_fallback = os.getenv("BOLTZ_ALLOW_CPU_FALLBACK", "0").strip().lower() in {
        "1",
        "true",
        "yes",
        "on",
    }
    summary = _log_gpu_runtime("predict")
    if summary.get("cuda_available"):
        return "gpu"
    if allow_cpu_fallback:
        logger.warning("[Boltz2] CUDA 不可用，按配置回退到 CPU")
        return "cpu"
    raise RuntimeError(
        f"Boltz2 GPU 不可用，拒绝回退到 CPU。 CUDA_VISIBLE_DEVICES={summary.get('cuda_visible_devices')!r}, torch={summary.get('torch_version')!r}"
    )


def _get_boltz_cache_dir() -> Path:
    default_cache_dir = Path(
        "/root/gpufree-data/amr_hunter_runtime/cache/boltz_model_cache"
    )
    return Path(os.getenv("BOLTZ_CACHE_DIR", str(default_cache_dir))).expanduser()


def _get_boltz_stable_job_cache_root() -> Path:
    configured = os.getenv("BOLTZ_STABLE_JOB_CACHE_DIR", "").strip()
    if configured:
        root = Path(configured).expanduser()
    else:
        root = Path("/root/gpufree-data/amr_hunter_runtime/cache/boltz_job_cache")
    root.mkdir(parents=True, exist_ok=True)
    return root


def _get_boltz_job_cache_root() -> Path:
    configured = os.getenv("BOLTZ_JOB_CACHE_DIR", "").strip()
    if configured:
        root = Path(configured).expanduser()
    else:
        root = _get_boltz_stable_job_cache_root()
    root.mkdir(parents=True, exist_ok=True)
    return root


def _job_key_from_payload(payload: Dict[str, Any]) -> str:
    boltz_input = _build_boltz_input(payload, attach_cached_msas=False)
    canonical = json.dumps(boltz_input, ensure_ascii=False, sort_keys=True)
    return hashlib.md5(canonical.encode("utf-8")).hexdigest()


def _job_dir_from_payload(payload: Dict[str, Any]) -> Path:
    job_key = _job_key_from_payload(payload)
    stable_dir = _get_boltz_stable_job_cache_root() / job_key
    if stable_dir.exists():
        return stable_dir
    active_dir = _get_boltz_job_cache_root() / job_key
    if active_dir.exists():
        return active_dir
    return stable_dir


def _prediction_output_root(input_file: Path, job_dir: Path) -> Path:
    return job_dir / f"boltz_results_{input_file.stem}"


def _persistent_runner_target() -> int:
    requested = max(
        1, int(os.getenv("BOLTZ_PERSISTENT_RUNNERS", str(_predict_concurrency)))
    )
    return max(1, min(requested, _predict_concurrency))


class _PersistentBoltzRunner:

    def __init__(self, runner_id: int) -> None:
        self.runner_id = runner_id
        self.cache_dir = _get_boltz_cache_dir()
        self.mol_dir = self.cache_dir / "mols"
        self.num_workers = max(0, int(os.getenv("BOLTZ_NUM_WORKERS", "0")))
        self.max_msa_seqs = int(os.getenv("BOLTZ_MAX_MSA_SEQS", "8192"))
        self.num_subsampled_msa = int(os.getenv("BOLTZ_NUM_SUBSAMPLED_MSA", "1024"))
        self.recycling_steps = int(os.getenv("BOLTZ_RECYCLING_STEPS", "3"))
        self.sampling_steps = int(os.getenv("BOLTZ_SAMPLING_STEPS", "200"))
        self.diffusion_samples = int(os.getenv("BOLTZ_DIFFUSION_SAMPLES", "1"))
        self.max_parallel_samples = int(os.getenv("BOLTZ_MAX_PARALLEL_SAMPLES", "5"))
        self.disable_kernels = _env_flag("BOLTZ_DISABLE_KERNELS", True)
        self.subsample_msa = _env_flag("BOLTZ_SUBSAMPLE_MSA", True)
        self.write_embeddings = _env_flag("BOLTZ_WRITE_EMBEDDINGS", False)
        self.output_format = (
            os.getenv("BOLTZ_OUTPUT_FORMAT", "pdb").strip().lower() or "pdb"
        )
        if self.output_format not in {"pdb", "mmcif"}:
            logger.warning("[Boltz2] 非法输出格式 %s，回退到 pdb", self.output_format)
            self.output_format = "pdb"
        self.accelerator = _resolve_accelerator()
        self.precision: Any = "bf16-mixed" if self.accelerator == "gpu" else 32
        self.runtime_dir = (
            _get_boltz_job_cache_root() / "__service_runtime" / f"runner_{runner_id}"
        )
        self.runtime_dir.mkdir(parents=True, exist_ok=True)
        self.model_module: Any = None
        self.trainer: Any = None
        self._initialize()

    def _initialize(self) -> None:
        import torch
        from pytorch_lightning import Trainer
        from boltz.data.write.writer import BoltzWriter
        from boltz.main import (
            Boltz2DiffusionParams,
            BoltzSteeringParams,
            MSAModuleArgs,
            PairformerArgsV2,
            download_boltz2,
        )
        from boltz.model.models.boltz2 import Boltz2

        started_at = time.perf_counter()
        download_boltz2(self.cache_dir)
        diffusion_params = Boltz2DiffusionParams()
        diffusion_params.step_scale = float(os.getenv("BOLTZ_STEP_SCALE", "1.5"))
        pairformer_args = PairformerArgsV2()
        msa_args = MSAModuleArgs(
            subsample_msa=self.subsample_msa,
            num_subsampled_msa=self.num_subsampled_msa,
            use_paired_feature=True,
        )
        steering_args = BoltzSteeringParams()
        predict_args = {
            "recycling_steps": self.recycling_steps,
            "sampling_steps": self.sampling_steps,
            "diffusion_samples": self.diffusion_samples,
            "max_parallel_samples": self.max_parallel_samples,
            "write_confidence_summary": True,
            "write_full_pae": False,
            "write_full_pde": False,
        }
        checkpoint_path = self.cache_dir / "boltz2_conf.ckpt"
        self.model_module = Boltz2.load_from_checkpoint(
            checkpoint_path,
            strict=True,
            predict_args=predict_args,
            map_location="cpu",
            diffusion_process_args=asdict(diffusion_params),
            ema=False,
            use_kernels=not self.disable_kernels,
            pairformer_args=asdict(pairformer_args),
            msa_args=asdict(msa_args),
            steering_args=asdict(steering_args),
        )
        self.model_module.eval()
        if hasattr(self.model_module, "freeze"):
            self.model_module.freeze()
        if self.accelerator == "gpu":
            self.model_module = self.model_module.to(torch.device("cuda"))
        dummy_writer = BoltzWriter(
            data_dir=str(self.runtime_dir / "processed" / "structures"),
            output_dir=str(self.runtime_dir / "predictions"),
            output_format=self.output_format,
            boltz2=True,
            write_embeddings=self.write_embeddings,
        )
        self.trainer = Trainer(
            default_root_dir=str(self.runtime_dir),
            strategy="auto",
            callbacks=[dummy_writer],
            accelerator=self.accelerator,
            devices=1,
            precision=self.precision,
            logger=False,
            enable_checkpointing=False,
            enable_progress_bar=False,
            enable_model_summary=False,
        )
        logger.info(
            "[Boltz2] 常驻运行时已就绪: runner=%d init=%.2fs accelerator=%s num_workers=%d",
            self.runner_id,
            time.perf_counter() - started_at,
            self.accelerator,
            self.num_workers,
        )

    def predict(self, prepared: Dict[str, Any]) -> None:
        from boltz.data.module.inferencev2 import Boltz2InferenceDataModule
        from boltz.data.write.writer import BoltzWriter

        processed = prepared["processed"]
        if not processed.manifest.records:
            return
        pred_writer = BoltzWriter(
            data_dir=str(processed.targets_dir),
            output_dir=str(prepared["output_root"] / "predictions"),
            output_format=self.output_format,
            boltz2=True,
            write_embeddings=self.write_embeddings,
        )
        if self.trainer.callbacks:
            self.trainer.callbacks[0] = pred_writer
        else:
            self.trainer.callbacks = [pred_writer]
        data_module = Boltz2InferenceDataModule(
            manifest=processed.manifest,
            target_dir=processed.targets_dir,
            msa_dir=processed.msa_dir,
            mol_dir=self.mol_dir,
            num_workers=self.num_workers,
            constraints_dir=processed.constraints_dir,
            template_dir=processed.template_dir,
            extra_mols_dir=processed.extra_mols_dir,
        )
        self.trainer.predict(
            self.model_module, datamodule=data_module, return_predictions=False
        )


def _ensure_persistent_runtime_pool() -> None:
    global _predict_concurrency, _predict_slots, _runtime_pool, _runtime_pool_ready, _runtime_init_error, _runtime_runner_count
    if not _persistent_runtime_enabled:
        return
    with _runtime_pool_guard:
        if _runtime_pool_ready:
            if _runtime_init_error:
                raise RuntimeError(_runtime_init_error)
            return
        pool: Queue[_PersistentBoltzRunner] = Queue()
        target_runners = _persistent_runner_target()
        init_errors: List[str] = []
        for runner_id in range(target_runners):
            try:
                pool.put(_PersistentBoltzRunner(runner_id))
            except Exception as exc:
                init_errors.append(str(exc))
                logger.exception("[Boltz2] 常驻运行时初始化失败: runner=%d", runner_id)
                break
        _runtime_runner_count = pool.qsize()
        _runtime_pool = pool if _runtime_runner_count else None
        _runtime_pool_ready = True
        if _runtime_runner_count == 0:
            _runtime_init_error = (
                init_errors[-1] if init_errors else "未能初始化任何常驻 Boltz 运行时"
            )
            raise RuntimeError(_runtime_init_error)
        if _runtime_runner_count != _predict_concurrency:
            logger.warning(
                "[Boltz2] 常驻运行时数量调整: requested=%d ready=%d",
                _predict_concurrency,
                _runtime_runner_count,
            )
            _predict_concurrency = _runtime_runner_count
            _predict_slots = threading.BoundedSemaphore(_predict_concurrency)
        if init_errors:
            logger.warning(
                "[Boltz2] 仅部分常驻运行时初始化成功: ready=%d requested=%d last_error=%s",
                _runtime_runner_count,
                target_runners,
                init_errors[-1],
            )


def _acquire_persistent_runner() -> _PersistentBoltzRunner:
    _ensure_persistent_runtime_pool()
    if _runtime_pool is None:
        raise RuntimeError(_runtime_init_error or "Boltz 常驻运行时池不可用")
    return _runtime_pool.get()


def _release_persistent_runner(runner: _PersistentBoltzRunner) -> None:
    if _runtime_pool is not None:
        _runtime_pool.put(runner)


def _load_prepared_boltz_job(output_root: Path) -> Dict[str, Any]:
    from boltz.main import BoltzProcessedInput, Manifest, filter_inputs_structure

    manifest = Manifest.load(output_root / "processed" / "manifest.json")
    filtered_manifest = filter_inputs_structure(
        manifest=manifest, outdir=output_root, override=False
    )
    processed_dir = output_root / "processed"
    processed = BoltzProcessedInput(
        manifest=filtered_manifest,
        targets_dir=processed_dir / "structures",
        msa_dir=processed_dir / "msa",
        constraints_dir=(
            processed_dir / "constraints"
            if (processed_dir / "constraints").exists()
            else None
        ),
        template_dir=(
            processed_dir / "templates"
            if (processed_dir / "templates").exists()
            else None
        ),
        extra_mols_dir=(
            processed_dir / "mols" if (processed_dir / "mols").exists() else None
        ),
    )
    return {
        "processed": processed,
        "total_records": len(manifest.records),
        "pending_records": len(filtered_manifest.records),
    }


def _prepare_boltz_job_for_persistent_runtime(
    payload: Dict[str, Any], input_file: Path, output_dir: Path
) -> Dict[str, Any]:
    from boltz.main import (
        BoltzProcessedInput,
        Manifest,
        check_inputs,
        filter_inputs_structure,
        process_inputs,
    )

    use_msa_server = _env_flag("BOLTZ_USE_MSA_SERVER", True)
    msa_server_url = (
        os.getenv("BOLTZ_MSA_SERVER_URL", "").strip() or "https://api.colabfold.com"
    )
    msa_pairing_strategy = (
        os.getenv("BOLTZ_MSA_PAIRING_STRATEGY", "").strip() or "greedy"
    )
    preprocessing_threads = max(1, int(os.getenv("BOLTZ_PREPROCESSING_THREADS", "1")))
    max_msa_seqs = int(os.getenv("BOLTZ_MAX_MSA_SEQS", "8192"))
    max_prepare_attempts = max(1, int(os.getenv("BOLTZ_PREPARE_PROCESS_ATTEMPTS", "3")))
    retry_delay_seconds = max(
        0.0, float(os.getenv("BOLTZ_PREPARE_PROCESS_RETRY_DELAY_SECONDS", "3"))
    )
    msa_server_username = os.getenv("BOLTZ_MSA_USERNAME")
    msa_server_password = os.getenv("BOLTZ_MSA_PASSWORD")
    api_key_header = os.getenv("BOLTZ_MSA_API_KEY_HEADER")
    api_key_value = os.getenv("MSA_API_KEY_VALUE")
    output_root = _prediction_output_root(input_file, output_dir)
    output_root.mkdir(parents=True, exist_ok=True)
    manifest_path = output_root / "processed" / "manifest.json"
    if manifest_path.exists():
        try:
            prepared = _load_prepared_boltz_job(output_root)
            if int(prepared["total_records"]) > 0:
                return {
                    "output_root": output_root,
                    "processed": prepared["processed"],
                    "reused_processed": True,
                    "local_msa_generated": 0,
                    "total_records": int(prepared["total_records"]),
                    "pending_records": int(prepared["pending_records"]),
                }
            logger.warning(
                "[Boltz2] 发现空 processed manifest，删除后重建: out=%s", output_root
            )
            shutil.rmtree(output_root / "processed", ignore_errors=True)
        except Exception as exc:
            logger.warning(
                "[Boltz2] 已存在 processed 结果不可复用，回退到重建: out=%s error=%s",
                output_root,
                exc,
            )
            shutil.rmtree(output_root / "processed", ignore_errors=True)
    reused_processed = False
    local_msa_generated = _materialize_local_msas_for_input(
        payload, input_file, output_dir
    )
    cache_dir = _get_boltz_cache_dir()
    job_key = _job_key_from_payload(payload)
    prepared = None
    remote_msa_required = False
    remote_msa_queue_wait = 0.0
    for attempt in range(1, max_prepare_attempts + 1):
        resolved_boltz_input = _build_boltz_input(payload)
        input_file.write_text(
            "# boltz input — auto-generated by AMR-Hunter\n"
            + _yaml_dump(resolved_boltz_input),
            encoding="utf-8",
        )
        remote_msa_required = use_msa_server and _boltz_input_needs_remote_msa(
            resolved_boltz_input
        )
        input_paths = check_inputs(input_file)
        with _acquire_remote_msa_slot(
            job_key, required=remote_msa_required
        ) as remote_msa_queue_wait:
            process_inputs(
                data=input_paths,
                out_dir=output_root,
                ccd_path=cache_dir / "ccd.pkl",
                mol_dir=cache_dir / "mols",
                msa_server_url=msa_server_url,
                msa_pairing_strategy=msa_pairing_strategy,
                max_msa_seqs=max_msa_seqs,
                use_msa_server=use_msa_server,
                msa_server_username=msa_server_username,
                msa_server_password=msa_server_password,
                api_key_header=api_key_header,
                api_key_value=api_key_value,
                boltz2=True,
                preprocessing_threads=preprocessing_threads,
            )
        _persist_stable_sequence_msas(
            resolved_boltz_input, output_root, input_file.stem
        )
        prepared = _load_prepared_boltz_job(output_root)
        if int(prepared["total_records"]) > 0:
            break
        if attempt >= max_prepare_attempts:
            break
        logger.warning(
            "[Boltz2] 输入预处理结果为空，准备重试: job=%s attempt=%d/%d remote_msa=%s out=%s",
            job_key[:12],
            attempt,
            max_prepare_attempts,
            remote_msa_required,
            output_root,
        )
        shutil.rmtree(output_root, ignore_errors=True)
        output_root.mkdir(parents=True, exist_ok=True)
        if retry_delay_seconds > 0:
            time.sleep(retry_delay_seconds * attempt)
    if prepared is None or int(prepared["total_records"]) <= 0:
        raise RuntimeError(
            f"Boltz 输入预处理未产生任何可预测记录。 manifest={output_root / 'processed' / 'manifest.json'}"
        )
    return {
        "output_root": output_root,
        "processed": prepared["processed"],
        "reused_processed": reused_processed,
        "local_msa_generated": local_msa_generated,
        "remote_msa_required": remote_msa_required,
        "remote_msa_queue_wait": remote_msa_queue_wait,
        "total_records": int(prepared["total_records"]),
        "pending_records": int(prepared["pending_records"]),
    }


def _prepare_boltz_inputs_with_slot(
    payload: Dict[str, Any], input_file: Path, output_dir: Path, *, error_prefix: str
) -> Dict[str, Any]:
    job_key = _job_key_from_payload(payload)
    duration_score = _payload_duration_score(payload)
    stdout_buffer = io.StringIO()
    stderr_buffer = io.StringIO()
    with _get_job_lock(job_key):
        with _acquire_prepare_slot(job_key, duration_score) as prepare_queue_wait:
            try:
                with (
                    contextlib.redirect_stdout(stdout_buffer),
                    contextlib.redirect_stderr(stderr_buffer),
                ):
                    prepare_started_at = time.perf_counter()
                    prepared = _prepare_boltz_job_for_persistent_runtime(
                        payload, input_file, output_dir
                    )
                    prepare_finished_at = time.perf_counter()
            except Exception as exc:
                raise BoltzPredictError(
                    f"{error_prefix}: {exc}",
                    stdout=stdout_buffer.getvalue(),
                    stderr=stderr_buffer.getvalue(),
                ) from exc
    return {
        "job_key": job_key,
        "prepared": prepared,
        "prepare_queue_wait": prepare_queue_wait,
        "prepare_started_at": prepare_started_at,
        "prepare_finished_at": prepare_finished_at,
        "stdout": stdout_buffer.getvalue(),
        "stderr": stderr_buffer.getvalue(),
    }


def _run_boltz_affinity_for_output(output_root: Path, accelerator: str) -> None:
    from pytorch_lightning import Trainer
    import torch
    from boltz.data.module.inferencev2 import Boltz2InferenceDataModule
    from boltz.data.write.writer import BoltzAffinityWriter
    from boltz.main import (
        Boltz2DiffusionParams,
        BoltzSteeringParams,
        MSAModuleArgs,
        Manifest,
        PairformerArgsV2,
        filter_inputs_affinity,
    )
    from boltz.model.models.boltz2 import Boltz2

    manifest_path = output_root / "processed" / "manifest.json"
    if not manifest_path.exists():
        return
    manifest = Manifest.load(manifest_path)
    if not any((record.affinity for record in manifest.records)):
        return
    manifest_filtered = filter_inputs_affinity(
        manifest=manifest, outdir=output_root, override=False
    )
    if not manifest_filtered.records:
        return
    processed_dir = output_root / "processed"
    runtime_dir = output_root / "__affinity_runtime"
    runtime_dir.mkdir(parents=True, exist_ok=True)
    diffusion_params = Boltz2DiffusionParams()
    pairformer_args = PairformerArgsV2()
    msa_args = MSAModuleArgs(
        subsample_msa=_env_flag("BOLTZ_SUBSAMPLE_MSA", True),
        num_subsampled_msa=int(os.getenv("BOLTZ_NUM_SUBSAMPLED_MSA", "1024")),
        use_paired_feature=True,
    )
    steering_args = BoltzSteeringParams()
    steering_args.fk_steering = False
    steering_args.physical_guidance_update = False
    steering_args.contact_guidance_update = False
    predict_affinity_args = {
        "recycling_steps": 5,
        "sampling_steps": int(os.getenv("BOLTZ_SAMPLING_STEPS_AFFINITY", "200")),
        "diffusion_samples": int(os.getenv("BOLTZ_DIFFUSION_SAMPLES_AFFINITY", "3")),
        "max_parallel_samples": 1,
        "write_confidence_summary": False,
        "write_full_pae": False,
        "write_full_pde": False,
    }
    affinity_checkpoint = _get_boltz_cache_dir() / "boltz2_aff.ckpt"
    model_module = Boltz2.load_from_checkpoint(
        affinity_checkpoint,
        strict=True,
        predict_args=predict_affinity_args,
        map_location="cpu",
        diffusion_process_args=asdict(diffusion_params),
        ema=False,
        pairformer_args=asdict(pairformer_args),
        msa_args=asdict(msa_args),
        steering_args=asdict(steering_args),
        affinity_mw_correction=_env_flag("BOLTZ_AFFINITY_MW_CORRECTION", False),
    )
    model_module.eval()
    if accelerator == "gpu":
        model_module = model_module.to(torch.device("cuda"))
    pred_writer = BoltzAffinityWriter(
        data_dir=str(processed_dir / "structures"),
        output_dir=str(output_root / "predictions"),
    )
    data_module = Boltz2InferenceDataModule(
        manifest=manifest_filtered,
        target_dir=output_root / "predictions",
        msa_dir=processed_dir / "msa",
        mol_dir=_get_boltz_cache_dir() / "mols",
        num_workers=max(0, int(os.getenv("BOLTZ_NUM_WORKERS", "0"))),
        constraints_dir=(
            processed_dir / "constraints"
            if (processed_dir / "constraints").exists()
            else None
        ),
        template_dir=(
            processed_dir / "templates"
            if (processed_dir / "templates").exists()
            else None
        ),
        extra_mols_dir=(
            processed_dir / "mols" if (processed_dir / "mols").exists() else None
        ),
        override_method="other",
        affinity=True,
    )
    trainer = Trainer(
        default_root_dir=str(runtime_dir),
        strategy="auto",
        callbacks=[pred_writer],
        accelerator=accelerator,
        devices=1,
        precision="bf16-mixed" if accelerator == "gpu" else 32,
        logger=False,
        enable_checkpointing=False,
        enable_progress_bar=False,
        enable_model_summary=False,
    )
    trainer.predict(model_module, datamodule=data_module, return_predictions=False)


def _run_boltz_predict_via_callback(
    payload: Dict[str, Any], input_file: Path, output_dir: Path
) -> Dict[str, Any]:
    from boltz.main import predict as boltz_predict

    use_msa_server = os.getenv("BOLTZ_USE_MSA_SERVER", "1").strip().lower() not in {
        "0",
        "false",
        "no",
        "off",
    }
    msa_server_url = os.getenv("BOLTZ_MSA_SERVER_URL", "").strip()
    msa_pairing_strategy = os.getenv("BOLTZ_MSA_PAIRING_STRATEGY", "").strip()
    disable_kernels = os.getenv("BOLTZ_DISABLE_KERNELS", "1").strip().lower() not in {
        "0",
        "false",
        "no",
        "off",
    }
    output_root = _prediction_output_root(input_file, output_dir)
    cache_dir = _get_boltz_cache_dir()
    job_key = _job_key_from_payload(payload)
    accelerator = _resolve_accelerator()
    stdout_buffer = io.StringIO()
    stderr_buffer = io.StringIO()
    had_predictions = any(output_root.rglob("confidence*.json")) or any(
        output_root.rglob("affinity*.json")
    )
    local_msa_generated = _materialize_local_msas_for_input(
        payload, input_file, output_dir
    )
    logger.info(
        "[Boltz2] 调度预测: job=%s accelerator=%s max_concurrent_predicts=%d reused=%s local_msa_generated=%d",
        job_key[:12],
        accelerator,
        _predict_concurrency,
        had_predictions,
        local_msa_generated,
    )
    predict_started_at = time.perf_counter()
    with _predict_slots:
        slot_acquired_at = time.perf_counter()
        logger.info(
            "[Boltz2] 进入预测槽位: job=%s queue_wait=%.2fs",
            job_key[:12],
            slot_acquired_at - predict_started_at,
        )
        with _get_job_lock(job_key):
            try:
                with (
                    contextlib.redirect_stdout(stdout_buffer),
                    contextlib.redirect_stderr(stderr_buffer),
                ):
                    callback_started_at = time.perf_counter()
                    boltz_predict.callback(
                        data=str(input_file),
                        out_dir=str(output_dir),
                        cache=str(cache_dir),
                        devices=1,
                        accelerator=accelerator,
                        output_format="pdb",
                        use_msa_server=use_msa_server,
                        msa_server_url=msa_server_url or "https://api.colabfold.com",
                        msa_pairing_strategy=msa_pairing_strategy or "greedy",
                        no_kernels=disable_kernels,
                        preprocessing_threads=int(
                            os.getenv("BOLTZ_PREPROCESSING_THREADS", "1")
                        ),
                        num_workers=int(os.getenv("BOLTZ_NUM_WORKERS", "0")),
                        max_msa_seqs=int(os.getenv("BOLTZ_MAX_MSA_SEQS", "8192")),
                        num_subsampled_msa=int(
                            os.getenv("BOLTZ_NUM_SUBSAMPLED_MSA", "1024")
                        ),
                    )
                    callback_finished_at = time.perf_counter()
            except Exception as exc:
                stdout_text = stdout_buffer.getvalue()
                stderr_text = stderr_buffer.getvalue()
                raise BoltzPredictError(
                    f"boltz predict Python callback 失败: {exc}",
                    stdout=stdout_text,
                    stderr=stderr_text,
                ) from exc
    stdout_text = stdout_buffer.getvalue()
    stderr_text = stderr_buffer.getvalue()
    if "Failed to process" in stdout_text and "Missing MSA" in stdout_text:
        raise BoltzPredictError(
            f"boltz predict 未生成有效输入：缺少 MSA，且未成功通过 MSA server 补全。\nstdout:\n{stdout_text[-2000:]}\nstderr:\n{stderr_text[-2000:]}",
            stdout=stdout_text,
            stderr=stderr_text,
        )
    logger.info(
        "[Boltz2] Python predict 完成: job=%s reused=%s output=%s slot_wait=%.2fs callback=%.2fs total=%.2fs",
        job_key[:12],
        had_predictions,
        output_root,
        slot_acquired_at - predict_started_at,
        callback_finished_at - callback_started_at,
        time.perf_counter() - predict_started_at,
    )
    _persist_stable_sequence_msas(
        _build_boltz_input(payload, attach_cached_msas=False),
        output_root,
        input_file.stem,
    )
    metrics = _parse_boltz_output(output_root)
    return {
        "raw": metrics,
        "stdout": stdout_text,
        "stderr": stderr_text,
        "job_dir": str(output_dir),
        "prediction_dir": str(output_root),
        "reused_prediction": had_predictions,
    }


def _run_boltz_predict(
    payload: Dict[str, Any], input_file: Path, output_dir: Path
) -> Dict[str, Any]:
    if not _persistent_runtime_enabled:
        return _run_boltz_predict_via_callback(payload, input_file, output_dir)
    output_root = _prediction_output_root(input_file, output_dir)
    job_key = _job_key_from_payload(payload)
    accelerator = _resolve_accelerator()
    had_predictions = any(output_root.rglob("confidence*.json")) or any(
        output_root.rglob("affinity*.json")
    )
    logger.info(
        "[Boltz2] 调度预测: job=%s accelerator=%s max_concurrent_predicts=%d mode=persistent reused=%s",
        job_key[:12],
        accelerator,
        _predict_concurrency,
        had_predictions,
    )
    predict_started_at = time.perf_counter()
    prepare_result = _prepare_boltz_inputs_with_slot(
        payload, input_file, output_dir, error_prefix="boltz 持久化运行时输入处理失败"
    )
    prepared = prepare_result["prepared"]
    prepare_queue_wait = float(prepare_result["prepare_queue_wait"])
    prepare_started_at = float(prepare_result["prepare_started_at"])
    prepare_finished_at = float(prepare_result["prepare_finished_at"])
    stdout_buffer = io.StringIO(prepare_result["stdout"])
    stderr_buffer = io.StringIO(prepare_result["stderr"])
    logger.info(
        "[Boltz2] 输入处理完成: job=%s queue_wait=%.2fs process=%.2fs records=%d pending=%d reused_processed=%s local_msa_generated=%d remote_msa=%s remote_wait=%.2fs",
        job_key[:12],
        prepare_queue_wait,
        prepare_finished_at - prepare_started_at,
        int(prepared["total_records"]),
        int(prepared["pending_records"]),
        bool(prepared["reused_processed"]),
        int(prepared.get("local_msa_generated") or 0),
        bool(prepared.get("remote_msa_required")),
        float(prepared.get("remote_msa_queue_wait") or 0.0),
    )
    if int(prepared["pending_records"]) == 0 and had_predictions:
        _run_boltz_affinity_for_output(output_root, accelerator)
        metrics = _parse_boltz_output(output_root)
        stdout_text = stdout_buffer.getvalue()
        stderr_text = stderr_buffer.getvalue()
        logger.info(
            "[Boltz2] 常驻推理完成: job=%s reused=%s output=%s process=%.2fs slot_wait=0.00s predict=0.00s total=%.2fs runner=<cache>",
            job_key[:12],
            had_predictions,
            output_root,
            prepare_finished_at - prepare_started_at,
            time.perf_counter() - predict_started_at,
        )
        return {
            "raw": metrics,
            "stdout": stdout_text,
            "stderr": stderr_text,
            "job_dir": str(output_dir),
            "prediction_dir": str(output_root),
            "reused_prediction": had_predictions,
        }
    runner: Optional[_PersistentBoltzRunner] = None
    slot_wait_started_at = time.perf_counter()
    try:
        runner = _acquire_persistent_runner()
        slot_acquired_at = time.perf_counter()
        logger.info(
            "[Boltz2] 进入预测槽位: job=%s queue_wait=%.2fs runner=%d",
            job_key[:12],
            slot_acquired_at - slot_wait_started_at,
            runner.runner_id,
        )
        with _get_job_lock(job_key):
            with (
                contextlib.redirect_stdout(stdout_buffer),
                contextlib.redirect_stderr(stderr_buffer),
            ):
                runner_started_at = time.perf_counter()
                runner.predict(prepared)
                _run_boltz_affinity_for_output(output_root, accelerator)
                runner_finished_at = time.perf_counter()
    except Exception as exc:
        stdout_text = stdout_buffer.getvalue()
        stderr_text = stderr_buffer.getvalue()
        raise BoltzPredictError(
            f"boltz 常驻运行时预测失败: {exc}", stdout=stdout_text, stderr=stderr_text
        ) from exc
    finally:
        if runner is not None:
            _release_persistent_runner(runner)
    stdout_text = stdout_buffer.getvalue()
    stderr_text = stderr_buffer.getvalue()
    if "Failed to process" in stdout_text and "Missing MSA" in stdout_text:
        raise BoltzPredictError(
            f"boltz predict 未生成有效输入：缺少 MSA，且未成功通过 MSA server 补全。\nstdout:\n{stdout_text[-2000:]}\nstderr:\n{stderr_text[-2000:]}",
            stdout=stdout_text,
            stderr=stderr_text,
        )
    logger.info(
        "[Boltz2] 常驻推理完成: job=%s reused=%s output=%s process=%.2fs slot_wait=%.2fs predict=%.2fs total=%.2fs runner=%d",
        job_key[:12],
        had_predictions,
        output_root,
        prepare_finished_at - prepare_started_at,
        slot_acquired_at - slot_wait_started_at,
        runner_finished_at - runner_started_at,
        time.perf_counter() - predict_started_at,
        runner.runner_id if runner is not None else -1,
    )
    metrics = _parse_boltz_output(output_root)
    return {
        "raw": metrics,
        "stdout": stdout_text,
        "stderr": stderr_text,
        "job_dir": str(output_dir),
        "prediction_dir": str(output_root),
        "reused_prediction": had_predictions,
    }


def _sanitize_name(value: Any) -> str:
    text = str(value or "unknown")
    safe = [char if char.isalnum() or char in ("-", "_", ".") else "_" for char in text]
    collapsed = "".join(safe).strip("._")
    return collapsed or "unknown"


def _tar_add_bytes(archive: tarfile.TarFile, arcname: str, payload: bytes) -> None:
    info = tarfile.TarInfo(name=arcname)
    info.size = len(payload)
    info.mtime = int(time.time())
    info.mode = 420
    archive.addfile(info, io.BytesIO(payload))


def _persist_boltz_artifacts(
    payload: Dict[str, Any],
    input_file: Path,
    output_dir: Path,
    *,
    raw: Optional[Dict[str, Any]] = None,
    stdout: str = "",
    stderr: str = "",
    error_message: Optional[str] = None,
) -> Optional[str]:
    meta = payload.get("_meta") if isinstance(payload.get("_meta"), dict) else {}
    root_text = str(
        meta.get("artifacts_dir") or os.getenv("BOLTZ_ARTIFACTS_DIR", "")
    ).strip()
    record_root_text = str(
        meta.get("artifact_record_root") or os.getenv("BOLTZ_ARTIFACT_RECORD_ROOT", "")
    ).strip()
    if not root_text:
        return None
    root_dir = Path(root_text)
    root_dir.mkdir(parents=True, exist_ok=True)
    mutation_id = _sanitize_name(meta.get("mutation_id"))
    gene_name = _sanitize_name(meta.get("gene_name"))
    seq_hash = _sanitize_name(meta.get("seq_hash"))[:12]
    scenario = _sanitize_name(meta.get("scenario") or payload.get("scenario"))
    stamp = datetime.utcnow().strftime("%Y%m%dT%H%M%SZ")
    artifact_name = (
        f"mut_{mutation_id}_{gene_name}_{scenario}_{seq_hash}_{stamp}.tar.gz"
    )
    artifact_path = root_dir / artifact_name
    compact_bundle_members = 0
    with tarfile.open(artifact_path, "w:gz") as archive:
        archive.add(input_file, arcname="input.yaml", recursive=False)
        if output_dir.exists():
            for path in sorted(output_dir.rglob("*")):
                if not path.is_file():
                    continue
                relative = path.relative_to(output_dir).as_posix()
                if relative.startswith(("msa/", "lightning_logs/", "processed/")):
                    continue
                if (
                    "/msa/" in relative
                    or "/lightning_logs/" in relative
                    or "/processed/" in relative
                ):
                    continue
                if path.suffix.lower() not in {
                    ".json",
                    ".pdb",
                    ".cif",
                    ".yaml",
                    ".txt",
                    ".log",
                }:
                    continue
                archive.add(path, arcname=f"output/{relative}")
                compact_bundle_members += 1
        metadata = {
            "meta": meta,
            "scenario": payload.get("scenario"),
            "persisted_at": stamp,
            "success": error_message is None,
            "artifact_mode": "bundle_tar_gz",
            "artifact_name": artifact_name,
            "local_prediction_dir": str(output_dir),
            "bundle_output_members": compact_bundle_members,
        }
        if record_root_text:
            metadata["share_archive_hint"] = (
                f"{Path(record_root_text).with_suffix('.tar')}#{artifact_name}"
            )
        _tar_add_bytes(
            archive,
            "metadata.json",
            json.dumps(metadata, ensure_ascii=False, indent=2).encode("utf-8"),
        )
        if raw is not None:
            _tar_add_bytes(
                archive,
                "raw_metrics.json",
                json.dumps(raw, ensure_ascii=False, indent=2).encode("utf-8"),
            )
        if stdout:
            _tar_add_bytes(archive, "stdout.log", stdout.encode("utf-8"))
        if stderr:
            _tar_add_bytes(archive, "stderr.log", stderr.encode("utf-8"))
        if error_message:
            _tar_add_bytes(archive, "error.txt", error_message.encode("utf-8"))
    if record_root_text:
        return f"{Path(record_root_text).with_suffix('.tar')}#{artifact_name}"
    return str(artifact_path)


def _persist_boltz_artifacts_safe(
    payload: Dict[str, Any],
    input_file: Path,
    output_dir: Path,
    *,
    raw: Optional[Dict[str, Any]] = None,
    stdout: str = "",
    stderr: str = "",
    error_message: Optional[str] = None,
) -> Optional[str]:
    try:
        return _persist_boltz_artifacts(
            payload,
            input_file,
            output_dir,
            raw=raw,
            stdout=stdout,
            stderr=stderr,
            error_message=error_message,
        )
    except Exception as exc:
        level = logging.WARNING if _is_disk_quota_error(exc) else logging.ERROR
        logger.log(level, "[Boltz2] 跳过共享盘 artifact 持久化: %s", exc)
        return None


def _parse_boltz_output(output_dir: Path) -> Dict[str, Any]:
    raw: Dict[str, Any] = {}
    confidence_found = False
    affinity_found = False
    for json_file in output_dir.rglob("confidence*.json"):
        try:
            data = json.loads(json_file.read_text(encoding="utf-8"))
            raw.update(data)
            confidence_found = True
            logger.info("[Boltz2] 解析置信度文件：%s", json_file)
            break
        except Exception as exc:
            logger.warning("[Boltz2] 解析 %s 失败：%s", json_file, exc)
    for json_file in output_dir.rglob("affinity*.json"):
        try:
            data = json.loads(json_file.read_text(encoding="utf-8"))
            raw.update(data)
            affinity_found = True
            logger.info("[Boltz2] 解析亲和力文件：%s", json_file)
            break
        except Exception as exc:
            logger.warning("[Boltz2] 解析 %s 失败：%s", json_file, exc)
    if not confidence_found and (not affinity_found):
        created_files = sorted(
            (str(path) for path in output_dir.rglob("*") if path.is_file())
        )
        tail = "\n".join(created_files[-20:]) if created_files else "<no files>"
        raise RuntimeError(
            f"Boltz 输出目录中未找到 confidence 或 affinity JSON。\noutput_dir={output_dir}\nfiles:\n{tail}"
        )
    return raw


def _cuda_available() -> bool:
    try:
        import torch

        return torch.cuda.is_available()
    except ImportError:
        return False


def _extract_metrics(raw: Dict[str, Any]) -> Dict[str, Optional[float]]:

    def _get(*keys: str) -> Optional[float]:
        for k in keys:
            v = raw.get(k)
            if isinstance(v, (int, float)) and (not v != v):
                return float(v)
            if isinstance(v, dict):
                mean_v = v.get("mean") or v.get("value")
                if isinstance(mean_v, (int, float)):
                    return float(mean_v)
        return None

    return {
        "iptm": _get("iptm", "interface_tm", "chain_iptm"),
        "ptm": _get("ptm"),
        "plddt": _get("complex_plddt", "plddt", "confidence"),
        "complex_energy": _get("complex_energy", "total_energy", "system_energy"),
        "binding_affinity": _get("binding_affinity", "affinity"),
        "affinity_pred_value": _get("affinity_pred_value"),
        "pair_energy": _get("pair_energy", "delta_g"),
        "stability_score": _get("stability", "folding_energy"),
        "clash_score": _get("clash_score", "clash", "steric_clash"),
    }


def build_app(model_dir: Path):
    try:
        from fastapi import FastAPI, HTTPException
        from pydantic import BaseModel
    except ImportError:
        logger.error(
            "缺少 fastapi/pydantic，请安装：pip install fastapi uvicorn[standard]"
        )
        sys.exit(1)
    app = FastAPI(title="Boltz2 Structure Prediction Server", version="1.0")
    _model_dir_ref = model_dir

    @app.on_event("startup")
    async def _startup():
        global _ready, _startup_error
        _prefetch_weights(_model_dir_ref)
        if _ready and _persistent_runtime_enabled:
            try:
                _ensure_persistent_runtime_pool()
            except Exception as exc:
                _ready = False
                _startup_error = f"Boltz 常驻运行时初始化失败：{exc}"
                logger.exception(_startup_error)

    @app.get("/health")
    async def health():
        if _startup_error:
            raise HTTPException(status_code=503, detail=_startup_error)
        if not _ready:
            raise HTTPException(status_code=503, detail="Boltz2 初始化中，请稍候")
        return {"status": "ok", "model": "boltz2"}

    @app.get("/ready")
    async def ready():
        return await health()

    class PredictResponse(BaseModel):
        iptm: Optional[float] = None
        ptm: Optional[float] = None
        plddt: Optional[float] = None
        complex_energy: Optional[float] = None
        binding_affinity: Optional[float] = None
        affinity_pred_value: Optional[float] = None
        pair_energy: Optional[float] = None
        stability_score: Optional[float] = None
        clash_score: Optional[float] = None
        artifact_dir: Optional[str] = None

    class PrefetchResponse(BaseModel):
        job_key: str
        reused_processed: bool
        local_msa_generated: int
        total_records: int
        pending_records: int
        prediction_cached: bool = False

    def _predict_handler(payload: Dict[str, Any]) -> PredictResponse:
        if not _ready:
            detail = _startup_error or "服务未就绪"
            raise HTTPException(status_code=503, detail=detail)
        work_dir = _job_dir_from_payload(payload)
        work_dir.mkdir(parents=True, exist_ok=True)
        request_started_at = time.perf_counter()
        dispatch_submitted_at = None
        meta = payload.get("_meta") if isinstance(payload.get("_meta"), dict) else {}
        if isinstance(meta, dict):
            dispatch_submitted_at = meta.get("dispatch_submitted_at")
        dispatch_lag_seconds: Optional[float] = None
        if dispatch_submitted_at:
            try:
                dispatch_lag_seconds = (
                    datetime.now(timezone.utc)
                    - datetime.fromisoformat(str(dispatch_submitted_at))
                ).total_seconds()
            except ValueError:
                dispatch_lag_seconds = None
        logger.info(
            "[Boltz2] 请求进入服务端: job=%s dispatch_lag=%s dur_score=%s",
            _job_key_from_payload(payload)[:12],
            (
                f"{dispatch_lag_seconds:.2f}s"
                if dispatch_lag_seconds is not None
                else "<unknown>"
            ),
            meta.get("estimated_duration_score") if isinstance(meta, dict) else None,
        )
        try:
            build_started_at = time.perf_counter()
            input_file = _write_boltz_input(payload, work_dir)
            output_dir = _prediction_output_root(input_file, work_dir)
            build_finished_at = time.perf_counter()
            logger.info(
                "[Boltz2] 输入准备完成: job=%s build=%.2fs",
                _job_key_from_payload(payload)[:12],
                build_finished_at - build_started_at,
            )
            run_result = _run_boltz_predict(payload, input_file, work_dir)
            raw = run_result["raw"]
            metrics = _extract_metrics(raw)
            persist_started_at = time.perf_counter()
            artifact_dir = _persist_boltz_artifacts_safe(
                payload,
                input_file,
                output_dir,
                raw=raw,
                stdout=run_result.get("stdout", ""),
                stderr=run_result.get("stderr", ""),
            )
            logger.info(
                "[Boltz2] 请求完成: job=%s total=%.2fs persist=%.2fs",
                _job_key_from_payload(payload)[:12],
                time.perf_counter() - request_started_at,
                time.perf_counter() - persist_started_at,
            )
            return PredictResponse(**metrics, artifact_dir=artifact_dir)
        except BoltzPredictError as exc:
            artifact_dir = _persist_boltz_artifacts_safe(
                payload,
                input_file,
                output_dir,
                stdout=exc.stdout,
                stderr=exc.stderr,
                error_message=str(exc),
            )
            detail = str(exc)
            if artifact_dir:
                detail = f"{detail} [artifacts: {artifact_dir}]"
            logger.exception("[Boltz2] 预测失败：%s", detail)
            raise HTTPException(status_code=500, detail=detail)
        except Exception as exc:
            detail = str(exc)
            if "input_file" in locals() and "output_dir" in locals():
                artifact_dir = _persist_boltz_artifacts_safe(
                    payload, input_file, output_dir, error_message=detail
                )
                if artifact_dir:
                    detail = f"{detail} [artifacts: {artifact_dir}]"
            logger.exception("[Boltz2] 预测失败：%s", detail)
            raise HTTPException(status_code=500, detail=detail)

    async def _predict(request: Request) -> PredictResponse:
        try:
            payload = await request.json()
        except Exception as exc:
            raise HTTPException(status_code=422, detail=f"JSON 解析失败：{exc}")
        return await asyncio.to_thread(_predict_handler, payload)

    def _prefetch_handler(payload: Dict[str, Any]) -> PrefetchResponse:
        if not _ready:
            detail = _startup_error or "服务未就绪"
            raise HTTPException(status_code=503, detail=detail)
        job_key = _job_key_from_payload(payload)
        work_dir = _job_dir_from_payload(payload)
        work_dir.mkdir(parents=True, exist_ok=True)
        request_started_at = time.perf_counter()
        meta = payload.get("_meta") if isinstance(payload.get("_meta"), dict) else {}
        logger.info(
            "[Boltz2] 预取请求进入服务端: job=%s dur_score=%s",
            job_key[:12],
            meta.get("estimated_duration_score") if isinstance(meta, dict) else None,
        )
        try:
            build_started_at = time.perf_counter()
            input_file = _write_boltz_input(payload, work_dir)
            output_dir = _prediction_output_root(input_file, work_dir)
            build_finished_at = time.perf_counter()
            logger.info(
                "[Boltz2] 预取输入准备完成: job=%s build=%.2fs",
                job_key[:12],
                build_finished_at - build_started_at,
            )
            prepare_result = _prepare_boltz_inputs_with_slot(
                payload, input_file, work_dir, error_prefix="boltz 预取输入处理失败"
            )
            prepared = prepare_result["prepared"]
            prediction_cached = any(output_dir.rglob("confidence*.json")) or any(
                output_dir.rglob("affinity*.json")
            )
            logger.info(
                "[Boltz2] 预取完成: job=%s queue_wait=%.2fs process=%.2fs pending=%d reused_processed=%s local_msa_generated=%d remote_msa=%s remote_wait=%.2fs total=%.2fs",
                job_key[:12],
                float(prepare_result["prepare_queue_wait"]),
                float(prepare_result["prepare_finished_at"])
                - float(prepare_result["prepare_started_at"]),
                int(prepared["pending_records"]),
                bool(prepared["reused_processed"]),
                int(prepared.get("local_msa_generated") or 0),
                bool(prepared.get("remote_msa_required")),
                float(prepared.get("remote_msa_queue_wait") or 0.0),
                time.perf_counter() - request_started_at,
            )
            return PrefetchResponse(
                job_key=job_key,
                reused_processed=bool(prepared["reused_processed"]),
                local_msa_generated=int(prepared.get("local_msa_generated") or 0),
                total_records=int(prepared["total_records"]),
                pending_records=int(prepared["pending_records"]),
                prediction_cached=prediction_cached,
            )
        except BoltzPredictError as exc:
            detail = str(exc)
            logger.exception("[Boltz2] 预取失败：%s", detail)
            raise HTTPException(status_code=500, detail=detail)
        except Exception as exc:
            detail = str(exc)
            logger.exception("[Boltz2] 预取失败：%s", detail)
            raise HTTPException(status_code=500, detail=detail)

    async def _prefetch(request: Request) -> PrefetchResponse:
        try:
            payload = await request.json()
        except Exception as exc:
            raise HTTPException(status_code=422, detail=f"JSON 解析失败：{exc}")
        return await asyncio.to_thread(_prefetch_handler, payload)

    app.post("/predict", response_model=PredictResponse)(_predict)
    app.post("/v1/predict", response_model=PredictResponse)(_predict)
    app.post("/structure", response_model=PredictResponse)(_predict)
    app.post("/v1/structure", response_model=PredictResponse)(_predict)
    app.post("/prefetch", response_model=PrefetchResponse)(_prefetch)
    app.post("/v1/prefetch", response_model=PrefetchResponse)(_prefetch)
    app.post("/prepare", response_model=PrefetchResponse)(_prefetch)
    app.post("/v1/prepare", response_model=PrefetchResponse)(_prefetch)
    return app


def main() -> None:
    parser = argparse.ArgumentParser(description="Boltz2 结构预测 HTTP 服务")
    parser.add_argument(
        "--model-dir",
        type=Path,
        default=Path(os.getenv("MODEL_DIR", "/root/gpufree-data/models/boltz2")),
        help="模型权重目录",
    )
    parser.add_argument(
        "--port", type=int, default=int(os.getenv("MODEL_PORT", "8001"))
    )
    parser.add_argument("--host", default="0.0.0.0")
    args = parser.parse_args()
    global _model_dir
    _model_dir = args.model_dir
    logger.info("[Boltz2 Server] 模型目录：%s", args.model_dir)
    logger.info("[Boltz2 Server] 监听：%s:%d", args.host, args.port)
    vram_limit = os.getenv("GPU_MEMORY_LIMIT_GB", "44")
    logger.info("[Boltz2 Server] 显存软限制：%s GB", vram_limit)
    logger.info(
        "[Boltz2 Server] Job cache：active=%s stable=%s",
        _get_boltz_job_cache_root(),
        _get_boltz_stable_job_cache_root(),
    )
    logger.info(
        "[Boltz2 Server] MSA cache：sequence=%s use_server=%s server=%s",
        _get_boltz_stable_msa_cache_root(),
        _env_flag("BOLTZ_USE_MSA_SERVER", True),
        os.getenv("BOLTZ_MSA_SERVER_URL", "").strip() or "https://api.colabfold.com",
    )
    logger.info(
        "[Boltz2 Server] Local MSA：enabled=%s require_local=%s bin=%s db=%s",
        _local_msa_enabled(),
        _local_msa_require_local(),
        _get_local_msa_search_bin(),
        _get_local_msa_db_dir() or "<unset>",
    )
    logger.info("[Boltz2 Server] 最大并发预测数：%d", _predict_concurrency)
    logger.info("[Boltz2 Server] 最大并发输入处理数：%d", _prepare_concurrency)
    logger.info("[Boltz2 Server] 最大并发远端 MSA 数：%d", _remote_msa_concurrency)
    logger.info(
        "[Boltz2 Server] 输入处理调度：long_threshold=%.1f long_cap=%d age_priority=%.0fs",
        _prepare_long_job_threshold,
        _prepare_long_job_cap,
        _prepare_age_priority_seconds,
    )
    logger.info(
        "[Boltz2 Server] 常驻运行时：%s",
        "enabled" if _persistent_runtime_enabled else "disabled",
    )
    if _persistent_runtime_enabled:
        logger.info(
            "[Boltz2 Server] 常驻 runner 目标数：%d", _persistent_runner_target()
        )
    logger.info(
        "[Boltz2 Server] DataLoader workers：%d",
        max(0, int(os.getenv("BOLTZ_NUM_WORKERS", "0"))),
    )
    _configure_torch_matmul_precision()
    _log_gpu_runtime("startup")
    try:
        import uvicorn
    except ImportError:
        logger.error("缺少 uvicorn，请安装：pip install uvicorn[standard]")
        sys.exit(1)
    app = build_app(args.model_dir)
    uvicorn.run(app, host=args.host, port=args.port, workers=1, log_level="info")


if __name__ == "__main__":
    main()
