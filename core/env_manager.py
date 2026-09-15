import os
import logging
from pathlib import Path
from typing import Optional
import shutil

logger = logging.getLogger(__name__)


def is_running_in_container() -> bool:
    return os.path.exists("/.dockerenv") or os.getenv("DOCKER_CONTAINER") == "true"


def get_platform_type() -> str:
    return os.getenv("AMR_PLATFORM_TYPE", "auto")


class EnvManager:
    LOCAL_DATA_VAR = "AMR_LOCAL_DATA"
    SHARE_PATH_VAR = "AMR_SHARE_PATH"
    PLATFORM_TYPE_VAR = "AMR_PLATFORM_TYPE"
    CONTAINER_LOCAL_ROOT = Path("/app/data")
    CONTAINER_SHARE_ROOT = Path("/app/share_data")
    CONTAINER_LOGS_ROOT = Path("/app/logs")
    HOST_LOCAL_ROOT = Path("/root/amr_hunter/data")
    HOST_SHARE_ROOT = Path("/root/gpufree-share/amr_hunter")
    HOST_LOGS_ROOT = Path("/root/amr_hunter/logs")
    CANDIDATE_SHARE_PATHS_CONTAINER = [
        Path("/app/share_data"),
        Path("/mnt/shared"),
        Path("/data/shared"),
    ]
    CANDIDATE_SHARE_PATHS_HOST = [
        Path("/root/gpufree-share/amr_hunter"),
        Path("/mnt/amr_share"),
        Path("/data/amr_share"),
        Path("./shared_data"),
    ]
    REQUIRED_DIRS = ["db", "results", "reference", "structures", "ligands"]

    def __init__(self):
        self.is_container = is_running_in_container()
        self.platform_type = get_platform_type()
        if self.is_container or self.platform_type == "container":
            self.local_root = Path(
                os.getenv(self.LOCAL_DATA_VAR, str(self.CONTAINER_LOCAL_ROOT))
            )
            self.share_root = self._resolve_share_path_container()
            self.logs_root = Path(
                os.getenv("AMR_LOG_DIR", str(self.CONTAINER_LOGS_ROOT))
            )
            self.env_name = "container"
        else:
            self.local_root = Path(
                os.getenv(self.LOCAL_DATA_VAR, str(self.HOST_LOCAL_ROOT))
            )
            self.share_root = self._resolve_share_path_host()
            self.logs_root = Path(os.getenv("AMR_LOG_DIR", str(self.HOST_LOGS_ROOT)))
            self.env_name = "host"
        logger.info(
            f"初始化环境管理器 (环境: {self.env_name}, 容器: {self.is_container})"
        )
        logger.info(f"  本地数据根: {self.local_root}")
        logger.info(f"  共享数据根: {self.share_root}")
        logger.debug(f"  日志根: {self.logs_root}")
        self._init_directories()

    def _resolve_share_path_container(self) -> Path:
        if env_share := os.getenv(self.SHARE_PATH_VAR):
            logger.info(f"📌 使用用户指定的共享路径: {env_share}")
            return Path(env_share)
        for candidate in self.CANDIDATE_SHARE_PATHS_CONTAINER:
            if candidate.exists():
                logger.info(f"✅ 探测到共享路径: {candidate}")
                return candidate
        logger.info(f"ℹ️  未检测到共享路径，使用默认值: {self.CONTAINER_SHARE_ROOT}")
        return self.CONTAINER_SHARE_ROOT

    def _resolve_share_path_host(self) -> Path:
        if env_share := os.getenv(self.SHARE_PATH_VAR):
            logger.info(f"📌 使用用户指定的共享路径: {env_share}")
            return Path(env_share)
        for candidate in self.CANDIDATE_SHARE_PATHS_HOST:
            if candidate.exists():
                logger.info(f"✅ 探测到共享路径: {candidate}")
                return candidate
        logger.info(f"ℹ️  未检测到共享路径，使用默认值: {self.HOST_SHARE_ROOT}")
        return self.HOST_SHARE_ROOT

    def _init_directories(self):
        self.local_root.mkdir(parents=True, exist_ok=True)
        for subdir in self.REQUIRED_DIRS:
            path = self.local_root / subdir
            path.mkdir(parents=True, exist_ok=True)
            logger.debug(f"目录就绪: {path}")
        self.logs_root.mkdir(parents=True, exist_ok=True)
        try:
            self.share_root.mkdir(parents=True, exist_ok=True)
            logger.debug(f"共享数据目录就绪: {self.share_root}")
        except OSError as e:
            logger.warning(
                f"无法创建共享目录 {self.share_root}: {e}（继续使用本地路径）"
            )

    def get_local_root(self) -> Path:
        return self.local_root

    def get_share_path(self) -> Path:
        return self.share_root

    def get_logs_path(self) -> Path:
        return self.logs_root

    def get_db_path(self) -> Path:
        return self.local_root / "db" / "amr_mutations.db"

    def get_results_dir(self) -> Path:
        return self.local_root / "results"

    def get_reference_dir(self) -> Path:
        return self.local_root / "reference"

    def get_structures_dir(self) -> Path:
        return self.local_root / "structures"

    def get_ligands_dir(self) -> Path:
        return self.local_root / "ligands"

    def sync_to_share(self, source_path: Path, subdir: str = "") -> bool:
        try:
            if not source_path.exists():
                logger.warning(f"源路径不存在: {source_path}")
                return False
            if subdir:
                target_root = self.share_root / subdir
            else:
                target_root = self.share_root
            target_root.mkdir(parents=True, exist_ok=True)
            if source_path.is_file():
                target_path = target_root / source_path.name
                shutil.copy2(source_path, target_path)
            else:
                target_path = target_root / source_path.name
                if target_path.exists():
                    shutil.rmtree(target_path)
                shutil.copytree(source_path, target_path)
            logger.info(f"已同步到共享存储: {source_path} → {target_path}")
            return True
        except Exception as e:
            logger.error(f"同步到共享存储失败: {e}")
            return False

    def restore_from_share(self, subdir: str, target_path: Path) -> bool:
        try:
            source_path = self.share_root / subdir
            if not source_path.exists():
                logger.debug(f"共享源不存在，跳过恢复: {source_path}")
                return False
            if target_path.exists():
                logger.debug(f"本地副本已存在，保留本地版本: {target_path}")
                return False
            if source_path.is_file():
                target_path.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source_path, target_path)
            else:
                shutil.copytree(source_path, target_path)
            logger.info(f"已从共享存储恢复: {source_path} → {target_path}")
            return True
        except Exception as e:
            logger.error(f"从共享存储恢复失败: {e}")
            return False

    def ensure_directory_exists(self, path: Path) -> bool:
        try:
            path.mkdir(parents=True, exist_ok=True)
            return True
        except Exception as e:
            logger.error(f"无法创建目录 {path}: {e}")
            return False

    def fallback_to_local(self, share_path: Path) -> Path:
        if not share_path.exists():
            logger.warning(f"共享路径不可用: {share_path}，降级使用本地路径")
            relative = (
                share_path.relative_to(self.share_root)
                if share_path.is_relative_to(self.share_root)
                else share_path.name
            )
            local_fallback = self.local_root / relative
            local_fallback.mkdir(parents=True, exist_ok=True)
            return local_fallback
        return share_path

    def get_environment_info(self) -> dict:
        return {
            "environment": self.env_name,
            "is_container": self.is_container,
            "platform_type": self.platform_type,
            "local_root": str(self.local_root),
            "share_root": str(self.share_root),
            "logs_root": str(self.logs_root),
            "db_path": str(self.get_db_path()),
            "results_dir": str(self.get_results_dir()),
            "reference_dir": str(self.get_reference_dir()),
            "structures_dir": str(self.get_structures_dir()),
            "ligands_dir": str(self.get_ligands_dir()),
        }

    def print_environment_info(self):
        info = self.get_environment_info()
        logger.info("环境信息:")
        for key, value in info.items():
            logger.info(f"  {key}: {value}")


_env_manager_instance: Optional[EnvManager] = None


def get_env_manager() -> EnvManager:
    global _env_manager_instance
    if _env_manager_instance is None:
        _env_manager_instance = EnvManager()
    return _env_manager_instance
