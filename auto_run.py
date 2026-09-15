import sys
import subprocess
import time
import logging
import os
import shutil
import tempfile
import re
from argparse import ArgumentParser
from pathlib import Path
from typing import Any, Dict, Optional
import json
from urllib.error import URLError
from urllib.request import Request, urlopen
import yaml

logging.basicConfig(
    level=logging.INFO, format="[%(asctime)s] [%(levelname)s] %(message)s"
)
logger = logging.getLogger(__name__)
LOCAL_DATA_DIR = Path(os.getenv("AMR_LOCAL_DATA", "/app/data"))
SHARE_DATA_DIR = Path(os.getenv("AMR_SHARE_PATH", "/app/share_data"))
LOG_DIR = Path(os.getenv("AMR_LOG_DIR", "/app/logs"))
DB_PATH = LOCAL_DATA_DIR / "db" / "amr_hunter.db"
RESULTS_DIR = LOCAL_DATA_DIR / "results"
REFERENCE_DIR = LOCAL_DATA_DIR / "reference"
STRUCTURES_DIR = LOCAL_DATA_DIR / "structures"
LIGANDS_DIR = LOCAL_DATA_DIR / "ligands"
INIT_SCRIPT_PATH = Path("/app/scripts/init_project.py")
DEFAULT_CONFIG_PATH = Path("/app/config/config.yaml")
PROJECT_PYTHON_BIN = os.getenv("AMR_PROJECT_PYTHON", sys.executable)
MAIN_SCRIPT_PATH = Path(os.getenv("AMR_MAIN_PATH", "/app/main.py"))
EVO2_URL = "http://evo2:8000"
BOLTZ2_URL = "http://boltz2:8001"


class HardwareRadar:

    @staticmethod
    def probe_nvidia_smi() -> dict:
        try:
            result = subprocess.run(
                [
                    "nvidia-smi",
                    "--query-gpu=index,name,memory.total,memory.free",
                    "--format=csv,noheader",
                ],
                capture_output=True,
                text=True,
                timeout=5,
            )
            if result.returncode != 0:
                raise RuntimeError(f"nvidia-smi failed: {result.stderr}")
            gpus = {}
            for line in result.stdout.strip().split("\n"):
                if not line.strip():
                    continue
                parts = [p.strip() for p in line.split(",")]
                gpu_id, name, total_mem, free_mem = (
                    parts[0],
                    parts[1],
                    parts[2],
                    parts[3],
                )
                gpus[gpu_id] = {
                    "name": name,
                    "total_memory": total_mem,
                    "free_memory": free_mem,
                }
            return {"status": "ok", "gpus": gpus}
        except Exception as e:
            logger.error(f"GPU probe failed: {e}")
            return {"status": "error", "error": str(e)}

    @staticmethod
    def check_docker_network() -> dict:
        checks = {}
        for service, url in [("evo2", EVO2_URL), ("boltz2", BOLTZ2_URL)]:
            try:
                result = subprocess.run(
                    ["curl", "-s", "-f", f"{url}/health"],
                    capture_output=True,
                    text=True,
                    timeout=3,
                )
                checks[service] = "ready" if result.returncode == 0 else "not ready"
            except Exception as e:
                checks[service] = f"error: {e}"
        return checks


class ServiceWarmup:
    MAX_RETRIES = 30
    RETRY_INTERVAL = 2

    @classmethod
    def warmup_all(cls):
        logger.info("🔥 开始服务预热...")
        logger.info(f"   等待 evo2 ({EVO2_URL})...")
        cls.warmup_service(EVO2_URL, "evo2")
        logger.info(f"   等待 boltz2 ({BOLTZ2_URL})...")
        cls.warmup_service(BOLTZ2_URL, "boltz2")
        logger.info("✅ 所有服务就绪")

    @classmethod
    def warmup_service(cls, url: str, service_name: str) -> bool:
        for attempt in range(1, cls.MAX_RETRIES + 1):
            try:
                result = subprocess.run(
                    ["curl", "-s", "-f", f"{url}/health"],
                    capture_output=True,
                    text=True,
                    timeout=3,
                )
                if result.returncode == 0:
                    logger.info(f"   ✅ {service_name} 就绪 (尝试 #{attempt})")
                    return True
            except Exception:
                pass
            if attempt < cls.MAX_RETRIES:
                logger.debug(
                    f"   ⏳ {service_name} 未就绪，{cls.RETRY_INTERVAL}s 后重试... (尝试 #{attempt})"
                )
                time.sleep(cls.RETRY_INTERVAL)
            else:
                logger.warning(
                    f"   ⚠️ {service_name} 在 {cls.MAX_RETRIES} 次尝试后仍未就绪"
                )
                return False
        return False


class DataSync:

    @staticmethod
    def sync_local_to_share(db_path: Optional[Path] = None):
        try:
            if not LOCAL_DATA_DIR.exists():
                logger.warning(f"本地数据目录不存在: {LOCAL_DATA_DIR}")
                return
            effective_db_path = db_path or DB_PATH
            SHARE_DATA_DIR.mkdir(parents=True, exist_ok=True)
            sync_items = [
                (effective_db_path, SHARE_DATA_DIR / "db"),
                (RESULTS_DIR, SHARE_DATA_DIR / "results"),
                (LOG_DIR, SHARE_DATA_DIR / "logs"),
            ]
            for src, dst in sync_items:
                if not src.exists():
                    logger.debug(f"跳过缺失源: {src}")
                    continue
                if src.is_file():
                    dst.mkdir(parents=True, exist_ok=True)
                    target = dst / src.name
                    shutil.copy2(src, target)
                    logger.info(f"✅ 同步文件: {src} → {target}")
                else:
                    if dst.exists():
                        shutil.rmtree(dst)
                    shutil.copytree(src, dst, dirs_exist_ok=True)
                    logger.info(f"✅ 同步目录: {src} → {dst}")
            if SHARE_DATA_DIR != LOCAL_DATA_DIR:
                logger.info(f"✅ 数据已同步到共享存储: {SHARE_DATA_DIR}")
        except Exception as e:
            logger.error(f"数据同步失败: {e}")

    @staticmethod
    def sync_share_to_local(db_path: Optional[Path] = None):
        try:
            if not SHARE_DATA_DIR.exists():
                logger.debug(f"共享数据目录不存在: {SHARE_DATA_DIR}")
                return
            effective_db_path = db_path or DB_PATH
            for item in ["results", "logs"]:
                src = SHARE_DATA_DIR / item
                dst = LOCAL_DATA_DIR / item
                if src.exists() and dst.exists():
                    logger.debug(f"保留本地版本，跳过覆盖: {dst}")
                elif src.exists():
                    shutil.copytree(src, dst, dirs_exist_ok=True)
                    logger.info(f"恢复: {src} → {dst}")
            shared_db = SHARE_DATA_DIR / "db" / effective_db_path.name
            if shared_db.exists() and (not effective_db_path.exists()):
                effective_db_path.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(shared_db, effective_db_path)
                logger.info(f"恢复数据库: {shared_db} → {effective_db_path}")
        except Exception as e:
            logger.warning(f"共享数据恢复失败: {e}")


class EvoModelInspector:
    MODEL_ENDPOINTS = ["/v1/models", "/models", "/model/info", "/health"]

    @staticmethod
    def _http_get_text(url: str, timeout_seconds: int = 3) -> Optional[str]:
        req = Request(url, headers={"User-Agent": "AMR-Hunter/1.0 (+model-check)"})
        try:
            with urlopen(req, timeout=max(1, timeout_seconds)) as response:
                return response.read().decode("utf-8", errors="ignore")
        except (URLError, OSError):
            return None

    @staticmethod
    def _detect_model_size(text: str) -> Optional[str]:
        payload = (text or "").lower()
        if "evo2_40b" in payload or re.search("\\b40b\\b", payload):
            return "40b"
        if "evo2_7b" in payload or re.search("\\b7b\\b", payload):
            return "7b"
        return None

    @classmethod
    def probe_loaded_model(cls) -> Dict[str, Optional[str]]:
        for endpoint in cls.MODEL_ENDPOINTS:
            url = f"{EVO2_URL}{endpoint}"
            text = cls._http_get_text(url)
            if not text:
                continue
            model_size = cls._detect_model_size(text)
            if model_size:
                return {"model_size": model_size, "endpoint": endpoint}
        return {"model_size": None, "endpoint": None}


def _append_suffix_to_path(path_text: str, suffix: str) -> str:
    path = Path(path_text)
    if path.suffix:
        return str(path.with_name(f"{path.stem}_{suffix}{path.suffix}"))
    return str(path.with_name(f"{path.name}_{suffix}"))


def _resolve_requested_model(cli_model: str) -> Optional[str]:
    if cli_model in {"7b", "40b"}:
        return cli_model
    env_model_name = str(os.getenv("EVO2_MODEL_NAME", "")).lower()
    if "40b" in env_model_name:
        return "40b"
    if "7b" in env_model_name:
        return "7b"
    env_model_size = str(os.getenv("EVO2_MODEL_SIZE", "")).lower()
    if env_model_size in {"7b", "40b"}:
        return env_model_size
    return None


def prepare_runtime_config(
    base_config_path: Path, requested_model: Optional[str]
) -> tuple[Path, Dict[str, Any], str]:
    with base_config_path.open("r", encoding="utf-8") as fp:
        config = yaml.safe_load(fp) or {}
    model_tag = requested_model or "unknown"
    suffix = f"evo{model_tag}"
    paths_cfg = config.setdefault("paths", {})
    pipeline_cfg = config.setdefault("pipeline", {})
    for key in ("local_db_path", "share_db_path"):
        current = paths_cfg.get(key)
        if current:
            paths_cfg[key] = _append_suffix_to_path(str(current), suffix)
    for key in ("export_csv_path", "log_path", "metrics_path"):
        current = pipeline_cfg.get(key)
        if current:
            pipeline_cfg[key] = _append_suffix_to_path(str(current), suffix)
    fd, temp_path = tempfile.mkstemp(prefix=f"amr_config_{suffix}_", suffix=".yaml")
    os.close(fd)
    runtime_config_path = Path(temp_path)
    with runtime_config_path.open("w", encoding="utf-8") as fp:
        yaml.safe_dump(config, fp, allow_unicode=True, sort_keys=False)
    return (runtime_config_path, config, model_tag)


class StageRunner:
    VALID_STAGES = ["evo", "boltz", "analyze", "all"]

    @staticmethod
    def run_init_project(config_path: Path) -> int:
        logger.info("📍 启动初始化: scripts/init_project.py")
        try:
            result = subprocess.run(
                [
                    PROJECT_PYTHON_BIN,
                    str(INIT_SCRIPT_PATH),
                    "--config",
                    str(config_path),
                ],
                timeout=None,
            )
            return result.returncode
        except Exception as e:
            logger.error(f"运行初始化失败: {e}")
            return 1

    @staticmethod
    def run_main_py(stage: str, config_path: Path) -> int:
        logger.info(f"📍 启动阶段: {stage}")
        try:
            result = subprocess.run(
                [
                    PROJECT_PYTHON_BIN,
                    str(MAIN_SCRIPT_PATH),
                    "--stage",
                    stage,
                    "--config",
                    str(config_path),
                ],
                timeout=None,
            )
            return result.returncode
        except Exception as e:
            logger.error(f"运行阶段失败 ({stage}): {e}")
            return 1

    @staticmethod
    def run_all_stages(config_path: Path) -> int:
        stages = ["evo", "boltz", "analyze"]
        for stage in stages:
            exit_code = StageRunner.run_main_py(stage, config_path=config_path)
            if exit_code != 0:
                logger.error(f"阶段 {stage} 失败 (exit code: {exit_code})")
                return exit_code
            logger.info(f"✅ 阶段 {stage} 完成")
        return 0


def parse_args() -> Any:
    parser = ArgumentParser(description="AMR-Hunter 容器指挥官")
    parser.add_argument(
        "--stage",
        choices=StageRunner.VALID_STAGES,
        default="all",
        help="执行阶段：evo|boltz|analyze|all",
    )
    parser.add_argument(
        "--evo-model",
        choices=["auto", "7b", "40b"],
        default="auto",
        help="期望的 Evo2 模型规格；auto 表示仅探测并尽力推断",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG_PATH,
        help="基础配置文件路径（将自动生成按模型隔离的运行时配置）",
    )
    return parser.parse_args()


def main():
    logger.info("=" * 70)
    logger.info("🚀 AMR-Hunter 容器指挥官启动")
    logger.info("=" * 70)
    runtime_config_path: Optional[Path] = None
    logger.info("📊 硬件探测...")
    radar = HardwareRadar()
    gpu_info = radar.probe_nvidia_smi()
    logger.info(f"   GPU 状态: {json.dumps(gpu_info, indent=4, ensure_ascii=False)}")
    network_checks = radar.check_docker_network()
    logger.info(
        f"   网络检查: {json.dumps(network_checks, indent=4, ensure_ascii=False)}"
    )
    args = parse_args()
    requested_model = _resolve_requested_model(args.evo_model)
    runtime_config_path, runtime_config, model_tag = prepare_runtime_config(
        base_config_path=args.config, requested_model=requested_model
    )
    runtime_db_path = Path(runtime_config["paths"]["local_db_path"])
    runtime_export_path = Path(runtime_config["pipeline"]["export_csv_path"])
    runtime_log_path = Path(runtime_config["pipeline"]["log_path"])
    logger.info(f"🧾 运行配置: {runtime_config_path}")
    logger.info(f"🧬 结果隔离标签: evo{model_tag}")
    logger.info("📁 初始化数据目录...")
    for path in [
        LOCAL_DATA_DIR,
        runtime_db_path.parent,
        RESULTS_DIR,
        REFERENCE_DIR,
        STRUCTURES_DIR,
        LIGANDS_DIR,
        LOG_DIR,
    ]:
        path.mkdir(parents=True, exist_ok=True)
        logger.debug(f"   ✓ {path}")
    logger.info("💾 恢复共享数据...")
    DataSync.sync_share_to_local(db_path=runtime_db_path)
    logger.info("🔥 服务预热...")
    ServiceWarmup.warmup_all()
    probe = EvoModelInspector.probe_loaded_model()
    detected_model = probe.get("model_size")
    if detected_model:
        logger.info(
            f"🧪 Evo2 当前模型探测: {detected_model} (endpoint={probe.get('endpoint')})"
        )
    else:
        logger.warning("⚠️ 未能从 Evo2 服务探测到 7b/40b 信息，将继续执行")
    stage = args.stage
    logger.info(f"📌 目标阶段: {stage}")
    if (
        args.evo_model in {"7b", "40b"}
        and detected_model
        and (detected_model != args.evo_model)
    ):
        logger.error(
            f"❌ 模型不匹配：命令指定 --evo-model {args.evo_model}，但探测到当前服务为 {detected_model}"
        )
        sys.exit(2)
    force_init = str(os.getenv("AMR_FORCE_INIT", "0")).lower() in {"1", "true", "yes"}
    should_init = (
        force_init
        or not runtime_db_path.exists()
        or runtime_db_path.stat().st_size == 0
    )
    if should_init:
        reason = "AMR_FORCE_INIT=1" if force_init else "数据库不存在或为空"
        logger.info(f"🧱 自动初始化触发: {reason}")
        init_code = StageRunner.run_init_project(config_path=runtime_config_path)
        if init_code != 0:
            logger.error(f"❌ 初始化失败 (exit code: {init_code})")
            sys.exit(init_code)
    else:
        logger.info(f"✅ 跳过初始化，检测到已有数据库: {runtime_db_path}")
    logger.info("⚙️  运行设定阶段...")
    if stage == "all":
        exit_code = StageRunner.run_all_stages(config_path=runtime_config_path)
    else:
        exit_code = StageRunner.run_main_py(stage, config_path=runtime_config_path)
    if exit_code != 0:
        logger.error(f"❌ 执行失败 (exit code: {exit_code})")
        sys.exit(exit_code)
    logger.info("📤 同步结果到共享存储...")
    DataSync.sync_local_to_share(db_path=runtime_db_path)
    logger.info("=" * 70)
    logger.info("✅ 容器指挥官执行完毕")
    logger.info("=" * 70)
    logger.info("")
    logger.info("📍 关键路径信息：")
    logger.info(f"   本地数据:  {LOCAL_DATA_DIR}")
    logger.info(f"   共享存储:  {SHARE_DATA_DIR}")
    logger.info(f"   结果文件:  {runtime_export_path}")
    logger.info(f"   数据库:    {runtime_db_path}")
    logger.info(f"   日志:      {runtime_log_path}")
    logger.info("")
    if runtime_config_path and runtime_config_path.exists():
        try:
            runtime_config_path.unlink(missing_ok=True)
        except OSError:
            pass
    sys.exit(0)


if __name__ == "__main__":
    main()
