from __future__ import annotations
import argparse
import logging
import os
import sys
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("download_models")
SHARE_MODELS = Path("/root/gpufree-share/amr_hunter_workspace/models")
EVO2_REPO = {"7b": "arcinstitute/evo2_7b", "40b": "arcinstitute/evo2_40b"}
BOLTZ2_REPO = "boltz-community/boltz-2"


def _require_huggingface_hub() -> None:
    try:
        import huggingface_hub
    except ImportError:
        logger.error(
            "缺少 huggingface_hub，请先安装：\n  pip install huggingface_hub\n  或：uv pip install huggingface_hub"
        )
        sys.exit(1)


def _login_if_needed(token: str | None) -> None:
    if not token:
        env_token = os.getenv("HF_TOKEN") or os.getenv("HUGGING_FACE_HUB_TOKEN")
        if env_token:
            token = env_token
    if token:
        from huggingface_hub import login

        login(token=token, add_to_git_credential=False)
        logger.info("[HF] 已使用 token 登录 Hugging Face Hub")
    else:
        logger.info("[HF] 未提供 token，以匿名模式下载（需要模型为公开仓库）")


def download_evo2(size: str, target_dir: Path, resume: bool = True) -> None:
    from huggingface_hub import snapshot_download

    repo_id = EVO2_REPO.get(size.lower())
    if not repo_id:
        raise ValueError(f"不支持的 Evo2 规格：{size}，可选：7b, 40b")
    dst = target_dir / f"evo2_{size}"
    dst.mkdir(parents=True, exist_ok=True)
    if (dst / "config.json").exists() and resume:
        logger.info("[Evo2] 权重已存在，跳过下载：%s", dst)
        return
    logger.info(
        "[Evo2] 开始下载 %s → %s（可能需要15–40分钟，取决于网速）", repo_id, dst
    )
    logger.info(
        "[Evo2] 文件大小参考：7B ≈ 14GB (fp16) | 40B ≈ 80GB (fp16) / 40GB (int8)"
    )
    snapshot_download(
        repo_id=repo_id,
        local_dir=str(dst),
        local_dir_use_symlinks=False,
        resume_download=True,
        ignore_patterns=["*.msgpack", "flax_model*"],
    )
    logger.info("[Evo2] 下载完成：%s", dst)


def download_boltz2(target_dir: Path, resume: bool = True) -> None:
    from huggingface_hub import snapshot_download

    dst = target_dir / "boltz2"
    dst.mkdir(parents=True, exist_ok=True)
    sentinels = ["boltz2.ckpt", "model.safetensors", "config.yaml", "config.json"]
    if any(((dst / s).exists() for s in sentinels)) and resume:
        logger.info("[Boltz2] 权重已存在，跳过下载：%s", dst)
        return
    logger.info("[Boltz2] 开始下载 %s → %s（约 3–6GB）", BOLTZ2_REPO, dst)
    try:
        snapshot_download(
            repo_id=BOLTZ2_REPO,
            local_dir=str(dst),
            local_dir_use_symlinks=False,
            resume_download=True,
        )
        logger.info("[Boltz2] 下载完成：%s", dst)
    except Exception as exc:
        logger.warning(
            "[Boltz2] HuggingFace 下载失败：%s\n  降级方案：boltz predict 首次运行时会自动下载权重到 ~/.boltz/\n  手动复制命令：cp -r ~/.boltz/boltz2 %s/",
            exc,
            dst,
        )


def main() -> int:
    parser = argparse.ArgumentParser(description="下载 AMR-Hunter 推理模型权重到共享盘")
    parser.add_argument(
        "--evo2-size",
        choices=["7b", "40b"],
        default="7b",
        help="Evo2 模型规格（default: 7b）",
    )
    parser.add_argument(
        "--only",
        choices=["evo2", "boltz2"],
        default=None,
        help="只下载指定模型（default: 下载两者）",
    )
    parser.add_argument(
        "--target-dir",
        type=Path,
        default=SHARE_MODELS,
        help=f"权重保存目录（default: {SHARE_MODELS}）",
    )
    parser.add_argument(
        "--hf-token",
        default=None,
        help="Hugging Face access token（受限模型需要，也可通过 HF_TOKEN 环境变量传入）",
    )
    parser.add_argument(
        "--no-resume", action="store_true", help="强制重新下载（忽略已存在的文件）"
    )
    args = parser.parse_args()
    _require_huggingface_hub()
    _login_if_needed(args.hf_token)
    resume = not args.no_resume
    target = args.target_dir
    target.mkdir(parents=True, exist_ok=True)
    logger.info("[存储] 权重目标目录：%s", target)
    try:
        if args.only in (None, "evo2"):
            download_evo2(args.evo2_size, target, resume=resume)
        if args.only in (None, "boltz2"):
            download_boltz2(target, resume=resume)
        logger.info("=" * 60)
        logger.info("✅ 模型下载完成")
        logger.info("   Evo2 %s → %s/evo2_%s/", args.evo2_size, target, args.evo2_size)
        logger.info("   Boltz2 → %s/boltz2/", target)
        logger.info("=" * 60)
        logger.info("下一步：运行 python run_pipeline.py 即可自动同步并启动服务")
        return 0
    except Exception as exc:
        logger.exception("下载失败：%s", exc)
        return 1


if __name__ == "__main__":
    sys.exit(main())
