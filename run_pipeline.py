from __future__ import annotations
import errno
import logging
import os
import re
import shlex
import shutil
import signal
import subprocess
import sys
import tarfile
import time
from argparse import ArgumentParser
from pathlib import Path
from typing import List, Optional
import yaml

SHARE_ROOT = Path("/root/gpufree-share/amr_hunter_workspace")
SHARE_CODE = SHARE_ROOT / "code"
SHARE_DATA = SHARE_ROOT / "amr_hunter_data"
SHARE_MODELS = SHARE_ROOT / "models"
SHARE_CACHE = SHARE_ROOT / "dependency_cache"
SHARE_BACKUP_ROOT = SHARE_ROOT / "backups"
FAST_ROOT = Path("/root/gpufree-data")
FAST_MODELS = FAST_ROOT / "models"
FAST_DB_DIR = FAST_ROOT / "db"
FAST_RUNTIME_ROOT = FAST_ROOT / "amr_hunter_runtime"
FAST_CACHE_DIR = FAST_RUNTIME_ROOT / "cache"
FAST_DEPENDENCY_CACHE = FAST_RUNTIME_ROOT / "dependency_cache"
FAST_TMP_DIR = FAST_RUNTIME_ROOT / "scratch" / "run_pipeline_tmp"
FAST_RESULTS = FAST_RUNTIME_ROOT / "results"
FAST_LOG_DIR = FAST_RUNTIME_ROOT / "logs"
SHARE_ROLLING_BACKUP_ARCHIVE = SHARE_BACKUP_ROOT / "amr_hunter_fast_state.tar.gz"
SHARE_DB = SHARE_DATA / "db" / "amr_hunter.db"
FAST_DB = FAST_DB_DIR / "amr_hunter.db"
CODE_DIR = Path(__file__).resolve().parent
CONFIG_PATH = CODE_DIR / "config" / "config.yaml"
MAIN_PY = CODE_DIR / "main.py"
PROJECT_VENV = Path(os.getenv("AMR_PROJECT_VENV", str(SHARE_CODE / ".venv")))
EVO2_VENV = Path(os.getenv("AMR_EVO2_VENV", str(SHARE_CODE / ".venv-evo2")))
BOLTZ2_VENV = Path(os.getenv("AMR_BOLTZ2_VENV", str(SHARE_CODE / ".venv-boltz2")))


def _resolve_python_bin(venv_dir: Path, override_env: str) -> str:
    override = os.getenv(override_env)
    if override:
        return override
    candidate = venv_dir / "bin" / "python"
    if candidate.exists():
        return str(candidate)
    return sys.executable


PROJECT_PYTHON_BIN = _resolve_python_bin(PROJECT_VENV, "AMR_PROJECT_PYTHON")
EVO2_PYTHON_BIN = _resolve_python_bin(EVO2_VENV, "AMR_EVO2_PYTHON")
BOLTZ2_PYTHON_BIN = _resolve_python_bin(BOLTZ2_VENV, "AMR_BOLTZ2_PYTHON")
EVO2_URL = "http://localhost:8000"
BOLTZ2_URL = "http://localhost:8001"
logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("run_pipeline")


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


_BOLTZ_TIMESTAMPED_CACHE_PATTERN = re.compile("^boltz_job_cache_.+_\\d{8}T\\d{6}$")


def _directory_size_bytes(path: Path) -> int:
    total = 0
    if not path.exists():
        return total
    for child in path.rglob("*"):
        if not child.is_file():
            continue
        try:
            total += child.stat().st_size
        except OSError:
            continue
    return total


def _directory_has_files(path: Path) -> bool:
    if not path.exists():
        return False
    for child in path.rglob("*"):
        if child.is_file():
            return True
    return False


def _resolve_boltz_job_cache_dir() -> Path:
    configured = str(os.getenv("BOLTZ_JOB_CACHE_DIR", "")).strip()
    if configured:
        return Path(configured).expanduser()
    return _resolve_boltz_stable_job_cache_dir()


def _resolve_boltz_stable_job_cache_dir() -> Path:
    configured = str(os.getenv("BOLTZ_STABLE_JOB_CACHE_DIR", "")).strip()
    if configured:
        return Path(configured).expanduser()
    return FAST_CACHE_DIR / "boltz_job_cache"


def _resolve_boltz_stable_msa_cache_dir() -> Path:
    configured = str(os.getenv("BOLTZ_STABLE_MSA_CACHE_DIR", "")).strip()
    if configured:
        return Path(configured).expanduser()
    return FAST_CACHE_DIR / "boltz_sequence_msas"


def _resolve_boltz_model_cache_dir() -> Path:
    configured = str(os.getenv("BOLTZ_CACHE_DIR", "")).strip()
    if configured:
        return Path(configured).expanduser()
    return FAST_CACHE_DIR / "boltz_model_cache"


def _prune_old_boltz_runtime_caches(active_cache_dir: Path) -> None:
    keep_recent = max(0, int(os.getenv("BOLTZ_JOB_CACHE_KEEP_RECENT_RUNS", "2")))
    cache_parent = active_cache_dir.parent
    if cache_parent != FAST_CACHE_DIR or not cache_parent.exists():
        return
    candidates = []
    for child in cache_parent.iterdir():
        if child == active_cache_dir or not child.is_dir():
            continue
        if not _BOLTZ_TIMESTAMPED_CACHE_PATTERN.match(child.name):
            continue
        try:
            candidates.append((child.stat().st_mtime, child))
        except OSError:
            continue
    candidates.sort(key=lambda item: (item[0], item[1].name), reverse=True)
    for _, stale_dir in candidates[keep_recent:]:
        reclaimed_gb = _directory_size_bytes(stale_dir) / 1024 / 1024 / 1024
        try:
            shutil.rmtree(stale_dir)
            logger.info(
                "[存储] 清理旧 Boltz 运行缓存: %s (%.2f GB)", stale_dir, reclaimed_gb
            )
        except Exception as exc:
            logger.warning("[存储] 清理旧 Boltz 运行缓存失败: %s (%s)", stale_dir, exc)


def _load_runtime_config(config_path: Path) -> dict:
    with config_path.open("r", encoding="utf-8") as fp:
        return yaml.safe_load(fp) or {}


def _apply_boltz_runtime_overrides(
    runtime_config: dict,
    *,
    service_predicts: Optional[int] = None,
    prewarm_runners: Optional[int] = None,
    remote_msa: Optional[int] = None,
    prefetch_workers: Optional[int] = None,
) -> dict:
    if not any(
        (
            value is not None
            for value in (
                service_predicts,
                prewarm_runners,
                remote_msa,
                prefetch_workers,
            )
        )
    ):
        return runtime_config
    pipeline_cfg = runtime_config.setdefault("pipeline", {})
    boltz_cfg = pipeline_cfg.setdefault("boltz", {})
    service_cfg = boltz_cfg.setdefault("service", {})
    dynamic_cfg = boltz_cfg.setdefault("dynamic_concurrency", {})
    msa_prefetch_cfg = service_cfg.setdefault("msa_prefetch", {})
    if service_predicts is not None:
        dynamic_cfg["service_max_concurrent_predicts"] = max(1, int(service_predicts))
    if prewarm_runners is not None:
        service_cfg["prewarm_runners"] = max(1, int(prewarm_runners))
    if remote_msa is not None:
        service_cfg["max_concurrent_remote_msa"] = max(1, int(remote_msa))
    if prefetch_workers is not None:
        msa_prefetch_cfg["max_workers"] = max(1, int(prefetch_workers))
    return runtime_config


def _write_runtime_config(config_path: Path, runtime_config: dict) -> None:
    with config_path.open("w", encoding="utf-8") as fp:
        yaml.safe_dump(runtime_config, fp, sort_keys=False, allow_unicode=True)


def _prepend_env_path(env: dict, key: str, path: Path) -> None:
    path_str = str(path)
    existing = env.get(key, "")
    parts = [part for part in existing.split(":") if part]
    if path_str in parts:
        parts.remove(path_str)
    env[key] = ":".join([path_str, *parts]) if parts else path_str


def _inject_cuda_library_paths(env: dict) -> None:
    candidates = (
        Path("/usr/local/cuda/lib64/stubs"),
        Path("/usr/local/cuda-11.8/compat"),
        Path("/usr/local/cuda-11.8/targets/x86_64-linux/lib/stubs"),
        Path("/usr/local/nvidia/lib64"),
    )
    for candidate in candidates:
        if not candidate.exists():
            continue
        if not (candidate / "libcuda.so").exists():
            continue
        _prepend_env_path(env, "LIBRARY_PATH", candidate)
        return


EVO2_GPU_MEMORY_LIMIT_GB = int(os.getenv("EVO2_GPU_MEMORY_LIMIT_GB", "44"))
BOLTZ2_GPU_MEMORY_LIMIT_GB = int(os.getenv("BOLTZ2_GPU_MEMORY_LIMIT_GB", "44"))


def _rsync_or_copy(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if shutil.which("rsync"):
        cmd = (
            ["rsync", "-a", "--progress", "--partial", str(src) + "/", str(dst) + "/"]
            if src.is_dir()
            else ["rsync", "-a", "--progress", "--partial", str(src), str(dst)]
        )
        result = subprocess.run(cmd, check=False)
        if result.returncode == 0:
            return
        logger.warning("rsync 返回非零码 %d，回退到 shutil", result.returncode)
    if src.is_dir():
        if dst.exists():
            shutil.rmtree(dst)
        shutil.copytree(src, dst)
    else:
        shutil.copy2(src, dst)


def _normalize_fast_backup_prefix(prefix: str) -> str:
    normalized = prefix.strip().strip("/")
    if normalized.startswith("fast/"):
        normalized = normalized[len("fast/") :]
    return normalized


def _restore_fast_state_from_share_backup(
    required_prefixes: Optional[List[str]] = None,
    dest_root: Path = FAST_ROOT,
    overwrite: bool = False,
) -> bool:
    archive_path = SHARE_ROLLING_BACKUP_ARCHIVE
    if not archive_path.exists():
        logger.warning(
            "[备份] 共享盘滚动备份不存在，无法恢复 FAST 状态: %s", archive_path
        )
        return False
    normalized_prefixes = None
    if required_prefixes:
        normalized_prefixes = [
            _normalize_fast_backup_prefix(prefix)
            for prefix in required_prefixes
            if _normalize_fast_backup_prefix(prefix)
        ]
        if not normalized_prefixes:
            return False
    restored_files = 0
    dest_root.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive_path, "r:gz") as archive:
        for member in archive.getmembers():
            if not member.isfile():
                continue
            member_name = member.name.strip("/")
            if not member_name.startswith("fast/"):
                continue
            rel_parts = Path(member_name).parts[1:]
            if not rel_parts or any((part in {"", ".", ".."} for part in rel_parts)):
                logger.warning("[备份] 跳过可疑备份成员: %s", member_name)
                continue
            rel_path = Path(*rel_parts)
            rel_path_str = rel_path.as_posix()
            if normalized_prefixes and (
                not any(
                    (
                        rel_path_str == prefix or rel_path_str.startswith(f"{prefix}/")
                        for prefix in normalized_prefixes
                    )
                )
            ):
                continue
            target_path = dest_root / rel_path
            if target_path.exists() and (not overwrite):
                continue
            source = archive.extractfile(member)
            if source is None:
                continue
            target_path.parent.mkdir(parents=True, exist_ok=True)
            with source, target_path.open("wb") as fp:
                shutil.copyfileobj(source, fp)
            restored_files += 1
    if restored_files:
        logger.info(
            "[备份] 已从共享盘滚动备份恢复 FAST 状态: %s → %s (%d files)",
            archive_path,
            dest_root,
            restored_files,
        )
        return True
    logger.info(
        "[备份] 滚动备份中没有可恢复的新文件: archive=%s prefixes=%s dest=%s",
        archive_path,
        normalized_prefixes or ["fast/"],
        dest_root,
    )
    return False


def _model_weights_present(model_dir: Path, model_name: str) -> bool:
    if model_name == "boltz2":
        sentinel_files = ["boltz2_aff.ckpt", "boltz2_conf.ckpt", "mols.tar"]
        return model_dir.exists() and all(
            ((model_dir / name).exists() for name in sentinel_files)
        )
    sentinel_files = [
        "config.json",
        "model.safetensors",
        "pytorch_model.bin",
        "model.pt",
    ]
    return model_dir.exists() and any(
        ((model_dir / name).exists() for name in sentinel_files)
    )


def sync_models_to_fast(evo2_size: str = "7b") -> None:
    FAST_MODELS.mkdir(parents=True, exist_ok=True)
    evo2_versioned = f"evo2_{evo2_size}"
    if (FAST_MODELS / evo2_versioned).exists() or (
        SHARE_MODELS / evo2_versioned
    ).exists():
        evo2_dir_name = evo2_versioned
    else:
        evo2_dir_name = "evo2"
    for model_name in (evo2_dir_name, "boltz2"):
        src = SHARE_MODELS / model_name
        dst = FAST_MODELS / model_name
        if _model_weights_present(dst, model_name):
            logger.info("[存储] %s 权重已在高速盘，跳过拷贝", model_name)
            continue
        if not src.exists():
            logger.warning(
                "[存储] %s 权重在高速盘缺失，且共享盘也不存在可恢复目录 %s",
                model_name,
                src,
            )
        else:
            logger.info(
                "[存储] 拷贝 %s 权重: %s → %s（可能需要数分钟）", model_name, src, dst
            )
            _rsync_or_copy(src, dst)
            logger.info("[存储] %s 权重拷贝完成", model_name)


def sync_db_to_fast() -> None:
    FAST_DB_DIR.mkdir(parents=True, exist_ok=True)
    if FAST_DB.exists():
        logger.info("[存储] 高速盘数据库已存在，按主副本保留: %s", FAST_DB)
        return
    if not SHARE_DB.exists():
        restored = _restore_fast_state_from_share_backup(
            required_prefixes=["db/amr_hunter.db", "db/amr_hunter.share_unsynced.db"],
            dest_root=FAST_ROOT,
            overwrite=False,
        )
        if restored and FAST_DB.exists():
            logger.info("[存储] 已从滚动备份恢复数据库到高速盘: %s", FAST_DB)
            return
        logger.warning("[存储] 共享盘数据库不存在: %s。将从头初始化。", SHARE_DB)
        return
    logger.info("[存储] 拷贝数据库到高速盘: %s → %s", SHARE_DB, FAST_DB)
    shutil.copy2(SHARE_DB, FAST_DB)
    logger.info("[存储] 数据库拷贝完成 (%.1f MB)", FAST_DB.stat().st_size / 1024 / 1024)


def sync_db_to_share() -> None:
    if not FAST_DB.exists():
        logger.warning("[存储] 高速盘数据库不存在，无法同步回共享盘")
        return
    SHARE_DB.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = SHARE_DB.with_suffix(".db.tmp")
    try:
        shutil.copy2(FAST_DB, tmp_path)
        tmp_path.replace(SHARE_DB)
        logger.info(
            "[存储] 数据库已同步回共享盘: %s (%.1f MB)",
            SHARE_DB,
            SHARE_DB.stat().st_size / 1024 / 1024,
        )
    except Exception as exc:
        if _is_disk_quota_error(exc):
            fallback_path = FAST_DB.with_name(
                f"{FAST_DB.stem}.share_unsynced{FAST_DB.suffix}"
            )
            shutil.copy2(FAST_DB, fallback_path)
            logger.warning(
                "[存储] 共享盘配额不足，跳过最终数据库共享同步；本地未同步快照保留在 %s",
                fallback_path,
            )
        else:
            logger.error("[存储] 同步数据库到共享盘失败: %s", exc)
            if tmp_path.exists():
                tmp_path.unlink(missing_ok=True)
            raise
        if tmp_path.exists():
            tmp_path.unlink(missing_ok=True)


def restore_runtime_state_to_fast_if_missing() -> None:
    missing_prefixes: List[str] = []
    if not _directory_has_files(FAST_RESULTS):
        missing_prefixes.append("results")
    if not _directory_has_files(FAST_LOG_DIR):
        missing_prefixes.append("logs")
    if not missing_prefixes:
        return
    _restore_fast_state_from_share_backup(
        required_prefixes=missing_prefixes, dest_root=FAST_ROOT, overwrite=False
    )


def sync_boltz_artifacts_to_share(local_root: Path, share_root: Path) -> None:
    if not local_root.exists():
        logger.info("[存储] 本地 Boltz artifacts 不存在，跳过归档: %s", local_root)
        return
    FAST_TMP_DIR.mkdir(parents=True, exist_ok=True)
    archive_root = FAST_RESULTS / "boltz_artifacts_archive"
    archive_root.mkdir(parents=True, exist_ok=True)
    share_archive = archive_root / f"{local_root.name}.tar"
    local_archive = FAST_TMP_DIR / f"{local_root.name}.tar"
    share_tmp = share_archive.with_name(f".{share_archive.name}.tmp")
    archive_ready = False
    logger.info(
        "[存储] 归档 Boltz artifacts 到 FAST: %s → %s", local_root, share_archive
    )
    try:
        local_archive.unlink(missing_ok=True)
        share_tmp.unlink(missing_ok=True)
        with tarfile.open(local_archive, "w") as archive:
            for path in sorted(local_root.rglob("*")):
                if not path.is_file():
                    continue
                arcname = f"{local_root.name}/{path.relative_to(local_root).as_posix()}"
                archive.add(path, arcname=arcname, recursive=False)
        shutil.copy2(local_archive, share_tmp)
        share_tmp.replace(share_archive)
        archive_ready = True
        logger.info(
            "[存储] Boltz artifacts 已归档到 FAST，本轮共享备份将统一收口: %s",
            share_archive,
        )
    finally:
        local_archive.unlink(missing_ok=True)
        if not archive_ready:
            share_tmp.unlink(missing_ok=True)
    try:
        shutil.rmtree(local_root)
        logger.info("[存储] 已清理本地 Boltz artifacts 目录: %s", local_root)
    except Exception as exc:
        logger.warning(
            "[存储] 清理本地 Boltz artifacts 目录失败: %s (%s)", local_root, exc
        )


def _add_path_to_archive(archive: tarfile.TarFile, src: Path, arc_root: str) -> None:
    if not src.exists():
        return
    if src.is_file():
        archive.add(src, arcname=arc_root, recursive=False)
        return
    for path in sorted(src.rglob("*")):
        if not path.is_file():
            continue
        archive.add(
            path,
            arcname=f"{arc_root}/{path.relative_to(src).as_posix()}",
            recursive=False,
        )


def update_share_rolling_backup() -> None:
    FAST_TMP_DIR.mkdir(parents=True, exist_ok=True)
    FAST_RESULTS.mkdir(parents=True, exist_ok=True)
    FAST_LOG_DIR.mkdir(parents=True, exist_ok=True)
    SHARE_BACKUP_ROOT.mkdir(parents=True, exist_ok=True)
    sources = [
        (FAST_DB, "fast/db/amr_hunter.db"),
        (
            FAST_DB_DIR / "amr_hunter.share_unsynced.db",
            "fast/db/amr_hunter.share_unsynced.db",
        ),
        (FAST_RESULTS, "fast/results"),
        (FAST_LOG_DIR, "fast/logs"),
        (SHARE_DATA / "reference", "share_inputs/reference"),
        (SHARE_DATA / "ligands", "share_inputs/ligands"),
        (SHARE_DATA / "structures", "share_inputs/structures"),
        (SHARE_CODE / "config", "code/config"),
    ]
    local_archive = (
        FAST_TMP_DIR / f"amr_hunter_fast_state_{time.strftime('%Y%m%dT%H%M%S')}.tar.gz"
    )
    share_tmp = SHARE_ROLLING_BACKUP_ARCHIVE.with_name(
        f".{SHARE_ROLLING_BACKUP_ARCHIVE.name}.tmp"
    )
    backup_ready = False
    logger.info("[备份] 生成共享盘滚动备份归档: %s", SHARE_ROLLING_BACKUP_ARCHIVE)
    try:
        local_archive.unlink(missing_ok=True)
        share_tmp.unlink(missing_ok=True)
        with tarfile.open(local_archive, "w:gz") as archive:
            for src, arc_root in sources:
                _add_path_to_archive(archive, src, arc_root)
        shutil.copy2(local_archive, share_tmp)
        share_tmp.replace(SHARE_ROLLING_BACKUP_ARCHIVE)
        backup_ready = True
        logger.info(
            "[备份] 共享盘滚动备份已更新: %s (%.1f MB)",
            SHARE_ROLLING_BACKUP_ARCHIVE,
            SHARE_ROLLING_BACKUP_ARCHIVE.stat().st_size / 1024 / 1024,
        )
    finally:
        local_archive.unlink(missing_ok=True)
        if not backup_ready:
            share_tmp.unlink(missing_ok=True)


def _build_service_env(
    gpu_memory_limit_gb: int, model_dir: Path, extra: dict | None = None
) -> dict:
    env = os.environ.copy()
    env["GPU_MEMORY_LIMIT_GB"] = str(gpu_memory_limit_gb)
    env["MODEL_DIR"] = str(model_dir)
    fast_hf_home = FAST_DEPENDENCY_CACHE / "huggingface"
    fast_xdg_cache = FAST_DEPENDENCY_CACHE / "xdg"
    fast_torch_home = FAST_DEPENDENCY_CACHE / "torch"
    fast_hf_home.mkdir(parents=True, exist_ok=True)
    fast_xdg_cache.mkdir(parents=True, exist_ok=True)
    fast_torch_home.mkdir(parents=True, exist_ok=True)
    env.setdefault("HF_HOME", str(fast_hf_home))
    env.setdefault("HF_HUB_CACHE", str(fast_hf_home / "hub"))
    env.setdefault("XDG_CACHE_HOME", str(fast_xdg_cache))
    env.setdefault("TORCH_HOME", str(fast_torch_home))
    env.setdefault("TMPDIR", str(FAST_TMP_DIR))
    env.setdefault("TMP", env["TMPDIR"])
    env.setdefault("TEMP", env["TMPDIR"])
    fraction = round(gpu_memory_limit_gb / 48.0, 2)
    env["PYTORCH_CUDA_ALLOC_CONF"] = f"max_split_size_mb:512"
    env["CUDA_VISIBLE_DEVICES"] = env.get("CUDA_VISIBLE_DEVICES", "0")
    FAST_TMP_DIR.mkdir(parents=True, exist_ok=True)
    _inject_cuda_library_paths(env)
    if extra:
        env.update({str(k): str(v) for k, v in extra.items()})
    return env


def start_service(
    name: str,
    start_cmd: str,
    model_dir: Path,
    gpu_memory_limit_gb: int,
    log_path: Path,
    extra_env: dict | None = None,
) -> subprocess.Popen:
    env = _build_service_env(gpu_memory_limit_gb, model_dir, extra_env)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_fh = log_path.open("a", buffering=1)
    args = shlex.split(start_cmd)
    logger.info("[服务] 启动 %s: %s", name, start_cmd)
    logger.info("[服务] 日志 → %s", log_path)
    proc = subprocess.Popen(
        args, stdout=log_fh, stderr=log_fh, env=env, preexec_fn=os.setsid
    )
    logger.info("[服务] %s PID=%d", name, proc.pid)
    return proc


def stop_service(proc: subprocess.Popen, name: str, timeout: int = 30) -> None:
    if proc.poll() is not None:
        logger.info("[服务] %s 已退出 (返回码=%d)", name, proc.returncode)
        return
    logger.info("[服务] 停止 %s (PID=%d)...", name, proc.pid)
    try:
        os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
    except ProcessLookupError:
        return
    try:
        proc.wait(timeout=timeout)
        logger.info("[服务] %s 已正常退出", name)
    except subprocess.TimeoutExpired:
        logger.warning("[服务] %s 未在 %ds 内退出，发送 SIGKILL", name, timeout)
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.wait()


def wait_for_service_healthy(
    url: str,
    name: str,
    proc: Optional[subprocess.Popen] = None,
    log_path: Optional[Path] = None,
    max_wait_seconds: int = 300,
    poll_interval: int = 5,
) -> bool:
    import socket
    import json
    import urllib.error
    import urllib.request
    from urllib.parse import urlparse

    parsed = urlparse(url)
    host = parsed.hostname or "localhost"
    port = parsed.port or 80
    health_paths = ("/health", "/v1/health", "/ready", "/")
    deadline = time.monotonic() + max_wait_seconds
    attempt = 0
    while time.monotonic() < deadline:
        attempt += 1
        if proc is not None and proc.poll() is not None:
            if log_path is not None:
                logger.error(
                    "[服务] %s 进程已提前退出 (返回码=%d)，请检查日志 %s",
                    name,
                    proc.returncode,
                    log_path,
                )
            else:
                logger.error(
                    "[服务] %s 进程已提前退出 (返回码=%d)", name, proc.returncode
                )
            return False
        try:
            with socket.create_connection((host, port), timeout=2):
                tcp_ok = True
        except OSError:
            tcp_ok = False
        if tcp_ok:
            for path in health_paths:
                try:
                    req = urllib.request.Request(f"{url}{path}", method="GET")
                    with urllib.request.urlopen(req, timeout=3) as resp:
                        if resp.status < 500:
                            logger.info(
                                "[服务] %s 已就绪（尝试 #%d，端点 %s）",
                                name,
                                attempt,
                                path,
                            )
                            return True
                except urllib.error.HTTPError as exc:
                    error_detail = ""
                    try:
                        payload = json.loads(exc.read().decode("utf-8"))
                        error_detail = str(payload.get("detail", "")).strip()
                    except Exception:
                        error_detail = ""
                    if (
                        exc.code == 503
                        and error_detail
                        and (error_detail != "模型加载中，请稍候")
                    ):
                        logger.error(
                            "[服务] %s 健康检查返回启动错误（端点 %s）: %s",
                            name,
                            path,
                            error_detail,
                        )
                        return False
                except Exception:
                    continue
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            break
        logger.info(
            "[服务] 等待 %s 就绪（尝试 #%d，剩余 %.0fs）...", name, attempt, remaining
        )
        time.sleep(min(poll_interval, remaining))
    logger.error("[服务] %s 在 %ds 内未能就绪", name, max_wait_seconds)
    return False


def run_stage(
    stage: str,
    config_path: Path,
    evo_model: str = "evo2_7b",
    limit: Optional[int] = None,
    gene_names: Optional[List[str]] = None,
    mutation_id_file: Optional[Path] = None,
    extra_env: Optional[dict] = None,
) -> int:
    cmd = [
        PROJECT_PYTHON_BIN,
        str(MAIN_PY),
        "--stage",
        stage,
        "--config",
        str(config_path),
        "--evo-model",
        evo_model,
    ]
    if limit is not None:
        cmd += ["--limit", str(limit)]
    for gene_name in gene_names or []:
        cmd += ["--gene", gene_name]
    if mutation_id_file is not None:
        cmd += ["--mutation-id-file", str(mutation_id_file)]
    logger.info("[流水线] 运行阶段 '%s': %s", stage, " ".join(cmd))
    env = os.environ.copy()
    if extra_env:
        env.update({k: str(v) for k, v in extra_env.items()})
    result = subprocess.run(cmd, check=False, env=env)
    if result.returncode != 0:
        logger.error("[流水线] 阶段 '%s' 失败，退出码=%d", stage, result.returncode)
    else:
        logger.info("[流水线] 阶段 '%s' 完成", stage)
    return result.returncode


def patch_config_for_fast_db(original_config_path: Path, fast_db_path: Path) -> Path:
    import tempfile
    import re

    content = original_config_path.read_text(encoding="utf-8")
    fast_export_path = FAST_RESULTS / "exports" / "amr_results.csv"
    fast_metrics_path = FAST_LOG_DIR / "batch_metrics.jsonl"
    fast_export_path.parent.mkdir(parents=True, exist_ok=True)
    fast_metrics_path.parent.mkdir(parents=True, exist_ok=True)
    content = re.sub("(local_db_path\\s*:\\s*).*", f"\\g<1>{fast_db_path}", content)
    content = re.sub("(share_db_path\\s*:\\s*).*", f"\\g<1>{fast_db_path}", content)
    content = re.sub(
        "(export_csv_path\\s*:\\s*).*", f"\\g<1>{fast_export_path}", content
    )
    content = re.sub("(metrics_path\\s*:\\s*).*", f"\\g<1>{fast_metrics_path}", content)
    tmp = tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".yaml",
        prefix="amr_config_fast_",
        delete=False,
        encoding="utf-8",
    )
    tmp.write(content)
    tmp.close()
    logger.info(
        "[配置] 临时配置已写入 %s（local_db_path/share_db_path → %s, export_csv_path → %s, metrics_path → %s）",
        tmp.name,
        fast_db_path,
        fast_export_path,
        fast_metrics_path,
    )
    return Path(tmp.name)


def request_shutdown(shutdown_command: str) -> bool:
    logger.warning("[关机] 流水线已完成，执行关机命令: %s", shutdown_command)
    try:
        result = subprocess.run(shlex.split(shutdown_command), check=False)
    except FileNotFoundError as exc:
        logger.error("[关机] 关机命令不存在: %s", exc)
        return False
    except Exception as exc:
        logger.error("[关机] 执行关机命令失败: %s", exc)
        return False
    if result.returncode != 0:
        logger.error("[关机] 关机命令返回非零退出码: %d", result.returncode)
        return False
    return True


def main() -> int:
    parser = ArgumentParser(
        description="AMR-Hunter 单容器调度脚本（L40S 顺序显存策略）"
    )
    parser.add_argument(
        "--stage",
        choices=["all", "evo", "boltz", "analyze"],
        default="all",
        help="要执行的流水线阶段（default: all）",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=CONFIG_PATH,
        help="config.yaml 路径（default: config/config.yaml）",
    )
    parser.add_argument(
        "--skip-model-sync",
        action="store_true",
        help="跳过模型权重同步（高速盘已有最新权重时使用）",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="限制每阶段处理的突变数量（测试运行用，如 --limit 10）",
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
    parser.add_argument(
        "--boltz-service-predicts",
        type=int,
        default=None,
        help="Temporarily override Boltz service_max_concurrent_predicts for this run.",
    )
    parser.add_argument(
        "--boltz-prewarm-runners",
        type=int,
        default=None,
        help="Temporarily override Boltz persistent prewarm_runners for this run.",
    )
    parser.add_argument(
        "--boltz-remote-msa",
        type=int,
        default=None,
        help="Temporarily override Boltz max_concurrent_remote_msa for this run.",
    )
    parser.add_argument(
        "--boltz-prefetch-workers",
        type=int,
        default=None,
        help="Temporarily override Boltz MSA prefetch worker count for this run.",
    )
    parser.add_argument(
        "--skip-db-sync", action="store_true", help="跳过数据库同步到高速盘（调试用）"
    )
    parser.add_argument(
        "--service-wait",
        type=int,
        default=300,
        help="等待服务就绪的最长秒数（default: 300）",
    )
    parser.add_argument(
        "--shutdown-on-complete",
        dest="shutdown_on_complete",
        action="store_true",
        default=None,
        help="流水线成功完成后执行关机命令（默认读取 config 中的 pipeline.shutdown_on_complete）",
    )
    parser.add_argument(
        "--no-shutdown-on-complete",
        dest="shutdown_on_complete",
        action="store_false",
        help="显式禁用完成后关机，覆盖 config 中的 pipeline.shutdown_on_complete",
    )
    parser.add_argument(
        "--shutdown-command",
        type=str,
        default=None,
        help="覆盖默认关机命令（例如 'shutdown -h now' 或 'systemctl poweroff'）",
    )
    parser.add_argument(
        "--evo2-size",
        choices=["7b", "20b", "40b"],
        default=os.getenv("EVO2_SIZE", "7b"),
        help="Evo2 模型规格（default: 7b）。20b 时显存 ~38GB（L40S 可用，但需要 Hopper 架构以正确支持 FP8；L40S 为 Ada Lovelace，精度可能回退）; 40b 时自动使用 int8 量化（~40GB VRAM）",
    )
    parser.add_argument(
        "--evo2-quantization",
        choices=["auto", "none", "int8", "int4"],
        default=os.getenv("EVO2_QUANTIZATION", "auto"),
        help="Evo2 量化模式（default: auto）。auto=目录名含 40b 自动选 int8；none=fp16 全精度；int8=精度损失 <1%%；int4=精度损失约 1-3%%，显存节省最多",
    )
    args = parser.parse_args()
    raw_config = _apply_boltz_runtime_overrides(
        _load_runtime_config(args.config),
        service_predicts=args.boltz_service_predicts,
        prewarm_runners=args.boltz_prewarm_runners,
        remote_msa=args.boltz_remote_msa,
        prefetch_workers=args.boltz_prefetch_workers,
    )
    pipeline_cfg = (
        raw_config.get("pipeline") or {} if isinstance(raw_config, dict) else {}
    )
    boltz_cfg = (
        pipeline_cfg.get("boltz") or {} if isinstance(pipeline_cfg, dict) else {}
    )
    boltz_service_cfg = (
        boltz_cfg.get("service") or {} if isinstance(boltz_cfg, dict) else {}
    )
    boltz_dynamic_cfg = (
        boltz_cfg.get("dynamic_concurrency") or {}
        if isinstance(boltz_cfg, dict)
        else {}
    )
    shutdown_on_complete = (
        args.shutdown_on_complete
        if args.shutdown_on_complete is not None
        else bool(pipeline_cfg.get("shutdown_on_complete", False))
    )
    shutdown_command = args.shutdown_command or str(
        pipeline_cfg.get("shutdown_command") or "shutdown -h now"
    )
    evo2_size = args.evo2_size
    evo_model_name = f"evo2_{evo2_size}"
    evo2_model_dir_env = os.getenv("EVO2_MODEL_DIR")
    boltz2_model_dir = Path(os.getenv("BOLTZ2_MODEL_DIR", str(FAST_MODELS / "boltz2")))
    evo2_quantization = args.evo2_quantization
    if evo2_quantization == "auto" and evo2_size == "40b":
        evo2_quantization = "int8"
        logger.info("[量化] 40B 模式，自动选择 int8（~40GB VRAM，精度损失 <1%%）")
        logger.info("[量化] 若需 int4（~22GB）请传入 --evo2-quantization int4")
    elif evo2_quantization == "auto" and evo2_size == "20b":
        evo2_quantization = "none"
        logger.info("[量化] 20B 模式，自动选择 none（原生 fp16，权重约 38GB）")
        logger.warning(
            "[量化/硬件] Evo2 20B 依赖 Transformer Engine FP8，需要 Hopper 架构 GPU（如 H100/H200）。L40S 为 Ada Lovelace，evo2 库将自动降级至 BF16/FP16，结果仍可用但可能有轻微差异。"
        )
    log_dir = FAST_LOG_DIR
    log_dir.mkdir(parents=True, exist_ok=True)
    evo2_proc: Optional[subprocess.Popen] = None
    boltz2_proc: Optional[subprocess.Popen] = None
    tmp_config: Optional[Path] = None
    run_id = time.strftime("run_%Y%m%dT%H%M%S")
    stable_boltz_job_cache_dir = _resolve_boltz_stable_job_cache_dir()
    configured_stable_msa_cache_dir = str(
        boltz_service_cfg.get("stable_msa_cache_dir") or ""
    ).strip()
    if configured_stable_msa_cache_dir:
        stable_boltz_msa_cache_dir = Path(configured_stable_msa_cache_dir).expanduser()
    else:
        stable_boltz_msa_cache_dir = _resolve_boltz_stable_msa_cache_dir()
    active_boltz_job_cache_dir = _resolve_boltz_job_cache_dir()
    boltz_model_cache_dir = _resolve_boltz_model_cache_dir()
    if active_boltz_job_cache_dir.parent == FAST_CACHE_DIR:
        FAST_CACHE_DIR.mkdir(parents=True, exist_ok=True)
        _prune_old_boltz_runtime_caches(active_boltz_job_cache_dir)
    stable_boltz_job_cache_dir.mkdir(parents=True, exist_ok=True)
    stable_boltz_msa_cache_dir.mkdir(parents=True, exist_ok=True)
    boltz_model_cache_dir.mkdir(parents=True, exist_ok=True)
    if active_boltz_job_cache_dir != stable_boltz_job_cache_dir:
        active_boltz_job_cache_dir.mkdir(parents=True, exist_ok=True)
    local_boltz_artifacts_root = FAST_RESULTS / "boltz_artifacts" / run_id
    share_boltz_artifacts_root = SHARE_DATA / "results" / "boltz_artifacts" / run_id
    boltz_artifacts_synced = False

    def _cleanup(sync_db: bool = True) -> None:
        nonlocal evo2_proc, boltz2_proc, tmp_config, boltz_artifacts_synced
        if evo2_proc:
            stop_service(evo2_proc, "evo2")
            evo2_proc = None
        if boltz2_proc:
            stop_service(boltz2_proc, "boltz2")
            boltz2_proc = None
        if not boltz_artifacts_synced:
            try:
                sync_boltz_artifacts_to_share(
                    local_boltz_artifacts_root, share_boltz_artifacts_root
                )
                boltz_artifacts_synced = True
            except Exception as exc:
                logger.warning("Boltz artifacts 批量同步失败: %s", exc)
        if sync_db:
            if str(
                os.getenv("AMR_ENABLE_LEGACY_SHARE_DB_SYNC", "")
            ).strip().lower() in {"1", "true", "yes", "on"}:
                try:
                    sync_db_to_share()
                except Exception as exc:
                    logger.error("数据库同步失败: %s", exc)
            else:
                logger.info(
                    "[存储] 已禁用共享盘散文件数据库同步；本轮仅更新滚动备份归档"
                )
            try:
                update_share_rolling_backup()
            except Exception as exc:
                logger.error("共享盘滚动备份失败: %s", exc)
        if tmp_config and tmp_config.exists():
            tmp_config.unlink(missing_ok=True)

    def _signal_handler(sig, frame):
        logger.warning("收到信号 %s，执行清理并退出...", sig)
        _cleanup(sync_db=True)
        sys.exit(1)

    signal.signal(signal.SIGTERM, _signal_handler)
    signal.signal(signal.SIGINT, _signal_handler)
    exit_code = 0
    should_request_shutdown = False
    try:
        logger.info("=" * 70)
        logger.info("AMR-Hunter 单容器部署流水线启动")
        logger.info("=" * 70)
        _check_gpu()
        if not args.skip_model_sync:
            sync_models_to_fast(evo2_size=evo2_size)
        else:
            logger.info("[存储] 跳过模型权重同步（--skip-model-sync）")
        if evo2_model_dir_env:
            evo2_model_dir = Path(evo2_model_dir_env)
        else:
            versioned_dir = FAST_MODELS / evo_model_name
            generic_dir = FAST_MODELS / "evo2"
            evo2_model_dir = versioned_dir if versioned_dir.exists() else generic_dir
        _services_dir = CODE_DIR / "services"
        evo2_start_cmd = os.getenv(
            "EVO2_START_CMD",
            f"{EVO2_PYTHON_BIN} {_services_dir / 'evo2_server.py'} --model-dir {evo2_model_dir} --port 8000 --quantization {evo2_quantization}",
        )
        boltz2_start_cmd = os.getenv(
            "BOLTZ2_START_CMD",
            f"{BOLTZ2_PYTHON_BIN} {_services_dir / 'boltz2_server.py'} --model-dir {boltz2_model_dir} --port 8001",
        )
        if not args.skip_db_sync:
            sync_db_to_fast()
        else:
            logger.info("[存储] 跳过数据库同步（--skip-db-sync）")
        restore_runtime_state_to_fast_if_missing()
        if FAST_DB.exists():
            tmp_config = patch_config_for_fast_db(args.config, FAST_DB)
            if any(
                (
                    value is not None
                    for value in (
                        args.boltz_service_predicts,
                        args.boltz_prewarm_runners,
                        args.boltz_remote_msa,
                        args.boltz_prefetch_workers,
                    )
                )
            ):
                patched_runtime_config = _apply_boltz_runtime_overrides(
                    _load_runtime_config(tmp_config),
                    service_predicts=args.boltz_service_predicts,
                    prewarm_runners=args.boltz_prewarm_runners,
                    remote_msa=args.boltz_remote_msa,
                    prefetch_workers=args.boltz_prefetch_workers,
                )
                _write_runtime_config(tmp_config, patched_runtime_config)
                logger.info(
                    "[配置] 已应用本轮 Boltz 并发 override 到临时配置: %s", tmp_config
                )
            active_config = tmp_config
        else:
            logger.warning("[配置] 高速盘 DB 不存在，使用原始配置（DB 路径不变）")
            active_config = args.config
        needs_evo = args.stage in ("all", "evo")
        needs_boltz = args.stage in ("all", "boltz")
        needs_analyze = args.stage in ("all", "analyze")
        if needs_evo:
            logger.info(
                "── Evo2 Stage ────────────────────────────────────────────────"
            )
            evo2_proc = start_service(
                name="evo2",
                start_cmd=evo2_start_cmd,
                model_dir=evo2_model_dir,
                gpu_memory_limit_gb=EVO2_GPU_MEMORY_LIMIT_GB,
                log_path=log_dir / "evo2_service.log",
                extra_env={
                    "MODEL_PORT": "8000",
                    "GPU_MEMORY_FRACTION": f"{EVO2_GPU_MEMORY_LIMIT_GB / 48:.2f}",
                },
            )
            if not wait_for_service_healthy(
                EVO2_URL,
                "evo2",
                proc=evo2_proc,
                log_path=log_dir / "evo2_service.log",
                max_wait_seconds=args.service_wait,
            ):
                logger.error("Evo2 服务未能在规定时间内就绪，中止流水线")
                _cleanup(sync_db=True)
                return 2
            rc = run_stage(
                "evo",
                active_config,
                evo_model=evo_model_name,
                limit=args.limit,
                gene_names=args.genes,
                mutation_id_file=args.mutation_id_file,
            )
            if rc != 0:
                exit_code = rc
            stop_service(evo2_proc, "evo2")
            evo2_proc = None
            logger.info("[显存] Evo2 已关闭，等待 GPU 显存释放...")
            time.sleep(5)
        if needs_boltz and exit_code == 0:
            logger.info(
                "── Boltz2 Stage ───────────────────────────────────────────────"
            )
            boltz2_proc = start_service(
                name="boltz2",
                start_cmd=boltz2_start_cmd,
                model_dir=boltz2_model_dir,
                gpu_memory_limit_gb=BOLTZ2_GPU_MEMORY_LIMIT_GB,
                log_path=log_dir / "boltz2_service.log",
                extra_env={
                    "MODEL_PORT": "8001",
                    "GPU_MEMORY_FRACTION": f"{BOLTZ2_GPU_MEMORY_LIMIT_GB / 48:.2f}",
                    "BOLTZ_CACHE_DIR": str(boltz_model_cache_dir),
                    "BOLTZ_ARTIFACTS_DIR": str(local_boltz_artifacts_root),
                    "BOLTZ_ARTIFACT_RECORD_ROOT": str(share_boltz_artifacts_root),
                    "BOLTZ_JOB_CACHE_DIR": str(active_boltz_job_cache_dir),
                    "BOLTZ_STABLE_JOB_CACHE_DIR": str(stable_boltz_job_cache_dir),
                    "BOLTZ_STABLE_MSA_CACHE_DIR": str(stable_boltz_msa_cache_dir),
                    "BOLTZ_SEQUENCE_MSA_CACHE": (
                        "1"
                        if bool(boltz_service_cfg.get("sequence_msa_cache", True))
                        else "0"
                    ),
                    "BOLTZ_USE_MSA_SERVER": (
                        "1"
                        if bool(boltz_service_cfg.get("use_msa_server", True))
                        else "0"
                    ),
                    "BOLTZ_MSA_SERVER_URL": str(
                        boltz_service_cfg.get("msa_server_url")
                        or "https://api.colabfold.com"
                    ),
                    "BOLTZ_MSA_PAIRING_STRATEGY": str(
                        boltz_service_cfg.get("msa_pairing_strategy") or "greedy"
                    ),
                    "BOLTZ_LOCAL_MSA_ENABLED": (
                        "1"
                        if bool(boltz_service_cfg.get("local_msa_enabled", False))
                        else "0"
                    ),
                    "BOLTZ_LOCAL_MSA_REQUIRE_LOCAL": (
                        "1"
                        if bool(boltz_service_cfg.get("local_msa_require_local", False))
                        else "0"
                    ),
                    "BOLTZ_LOCAL_MSA_BIN": str(
                        boltz_service_cfg.get("local_msa_bin") or "colabfold_search"
                    ),
                    "BOLTZ_LOCAL_MSA_DB_DIR": str(
                        boltz_service_cfg.get("local_msa_db_dir") or ""
                    ),
                    "BOLTZ_LOCAL_MSA_MMSEQS_BIN": str(
                        boltz_service_cfg.get("local_msa_mmseqs_bin") or ""
                    ),
                    "BOLTZ_LOCAL_MSA_EXTRA_ARGS": str(
                        boltz_service_cfg.get("local_msa_extra_args") or ""
                    ),
                    "BOLTZ_LOCAL_MSA_WORK_DIR": str(
                        boltz_service_cfg.get("local_msa_work_dir") or ""
                    ),
                    "BOLTZ_MAX_CONCURRENT_PREDICTS": str(
                        max(
                            1,
                            int(
                                boltz_dynamic_cfg.get(
                                    "service_max_concurrent_predicts",
                                    boltz_dynamic_cfg.get(
                                        "max_workers", boltz_cfg.get("workers", 1)
                                    ),
                                )
                            ),
                        )
                    ),
                    "BOLTZ_PERSISTENT_RUNTIME": (
                        "1"
                        if bool(boltz_service_cfg.get("persistent_runtime", True))
                        else "0"
                    ),
                    "BOLTZ_PERSISTENT_RUNNERS": str(
                        max(
                            1,
                            int(
                                boltz_service_cfg.get(
                                    "prewarm_runners",
                                    boltz_dynamic_cfg.get(
                                        "service_max_concurrent_predicts",
                                        boltz_dynamic_cfg.get(
                                            "max_workers", boltz_cfg.get("workers", 1)
                                        ),
                                    ),
                                )
                            ),
                        )
                    ),
                    "BOLTZ_NUM_WORKERS": str(
                        max(0, int(boltz_service_cfg.get("num_workers", 0)))
                    ),
                    "BOLTZ_PREPROCESSING_THREADS": str(
                        max(1, int(boltz_service_cfg.get("preprocessing_threads", 1)))
                    ),
                    "BOLTZ_MAX_MSA_SEQS": str(
                        max(1, int(boltz_service_cfg.get("max_msa_seqs", 8192)))
                    ),
                    "BOLTZ_NUM_SUBSAMPLED_MSA": str(
                        max(1, int(boltz_service_cfg.get("num_subsampled_msa", 1024)))
                    ),
                    "BOLTZ_SUBSAMPLE_MSA": (
                        "1"
                        if bool(boltz_service_cfg.get("subsample_msa", True))
                        else "0"
                    ),
                    "BOLTZ_MAX_CONCURRENT_PREPARES": str(
                        max(1, int(boltz_service_cfg.get("max_concurrent_prepares", 1)))
                    ),
                    "BOLTZ_MAX_CONCURRENT_REMOTE_MSA": str(
                        max(
                            1,
                            int(boltz_service_cfg.get("max_concurrent_remote_msa", 1)),
                        )
                    ),
                    "BOLTZ_PREPARE_PROCESS_ATTEMPTS": str(
                        max(
                            1, int(boltz_service_cfg.get("prepare_process_attempts", 3))
                        )
                    ),
                    "BOLTZ_PREPARE_PROCESS_RETRY_DELAY_SECONDS": str(
                        max(
                            0.0,
                            float(
                                boltz_service_cfg.get(
                                    "prepare_process_retry_delay_seconds", 3
                                )
                            ),
                        )
                    ),
                },
            )
            if not wait_for_service_healthy(
                BOLTZ2_URL,
                "boltz2",
                proc=boltz2_proc,
                log_path=log_dir / "boltz2_service.log",
                max_wait_seconds=args.service_wait,
            ):
                logger.error("Boltz2 服务未能在规定时间内就绪，中止流水线")
                _cleanup(sync_db=True)
                return 2
            rc = run_stage(
                "boltz",
                active_config,
                evo_model=evo_model_name,
                limit=args.limit,
                gene_names=args.genes,
                mutation_id_file=args.mutation_id_file,
                extra_env={
                    "BOLTZ_ARTIFACTS_DIR": str(local_boltz_artifacts_root),
                    "BOLTZ_ARTIFACT_RECORD_ROOT": str(share_boltz_artifacts_root),
                    "BOLTZ_JOB_CACHE_DIR": str(active_boltz_job_cache_dir),
                    "BOLTZ_STABLE_JOB_CACHE_DIR": str(stable_boltz_job_cache_dir),
                    "BOLTZ_STABLE_MSA_CACHE_DIR": str(stable_boltz_msa_cache_dir),
                },
            )
            if rc != 0:
                exit_code = rc
            stop_service(boltz2_proc, "boltz2")
            boltz2_proc = None
            try:
                sync_boltz_artifacts_to_share(
                    local_boltz_artifacts_root, share_boltz_artifacts_root
                )
                boltz_artifacts_synced = True
            except Exception as exc:
                logger.warning("Boltz artifacts 批量同步失败: %s", exc)
        if needs_analyze and exit_code == 0:
            logger.info(
                "── Analyze Stage ──────────────────────────────────────────────"
            )
            rc = run_stage("analyze", active_config, evo_model=evo_model_name)
            if rc != 0:
                exit_code = rc
        logger.info("=" * 70)
        logger.info("流水线阶段完成，同步数据到共享盘...")
        sync_db_to_share()
        if exit_code == 0:
            logger.info("✅ AMR-Hunter 流水线全部完成")
            if shutdown_on_complete:
                should_request_shutdown = True
                logger.warning("[关机] 已启用完成后关机，将在清理完成后执行")
        else:
            logger.warning("⚠️  流水线完成但有阶段失败（exit_code=%d）", exit_code)
    except Exception as exc:
        logger.exception("流水线发生未预期错误: %s", exc)
        exit_code = 1
    finally:
        _cleanup(sync_db=True)
    if should_request_shutdown:
        request_shutdown(shutdown_command)
    return exit_code


def _check_gpu() -> None:
    try:
        result = subprocess.run(
            [
                "nvidia-smi",
                "--query-gpu=index,name,memory.total,memory.free,utilization.gpu",
                "--format=csv,noheader",
            ],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
        if result.returncode == 0:
            for line in result.stdout.strip().splitlines():
                logger.info("[GPU] %s", line.strip())
        else:
            logger.warning("[GPU] nvidia-smi 失败: %s", result.stderr.strip())
    except FileNotFoundError:
        logger.warning("[GPU] nvidia-smi 未找到，跳过 GPU 检查")
    except OSError as exc:
        logger.warning("[GPU] nvidia-smi 不可执行，跳过 GPU 检查: %s", exc)
    except subprocess.TimeoutExpired:
        logger.warning("[GPU] nvidia-smi 超时")


if __name__ == "__main__":
    sys.exit(main())
