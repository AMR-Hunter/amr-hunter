from __future__ import annotations
import argparse
import importlib
import logging
import math
import os
import shutil
import site
import sys
import time
from pathlib import Path
from typing import Any, Dict, Optional

try:
    from pydantic import BaseModel
except ImportError:
    BaseModel = object
logger = logging.getLogger("evo2_server")
logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
_model = None
_tokenizer = None
_model_name: str = "evo2_7b"
_device: str = "cuda"
_quantization: str = "none"


class LikelihoodRequest(BaseModel):
    sequence: str


class LikelihoodResponse(BaseModel):
    log_likelihood: float
    model: str = _model_name


def _detect_quantization(model_dir: Path, requested: str) -> str:
    if requested != "auto":
        return requested
    dir_name = model_dir.name.lower()
    if "40b" in dir_name:
        logger.info("[量化] 检测到 40B 模型目录，自动启用 int8")
        return "int8"
    if "20b" in dir_name:
        logger.info("[量化] 检测到 20B 模型目录，使用 none（fp16，权重 ~38GB）")
        logger.warning(
            "[量化/硬件] Evo2 20B 依赖 Transformer Engine FP8 Linear，需要 Hopper 架构 GPU（H100/H200/GH200）。当前硬件若为 Ada Lovelace（L40S/RTX 4090）则 TE 将自动降级至 BF16/FP16，结果仍可用但可能存在轻微数值差异。"
        )
        return "none"
    config_file = model_dir / "config.json"
    if config_file.exists():
        try:
            import json

            cfg = json.loads(config_file.read_text(encoding="utf-8"))
            num_params = cfg.get("num_parameters", 0) or 0
            if num_params > 20000000000:
                logger.info(
                    "[量化] 参数量 %.0fB 超过 20B，自动启用 int8",
                    num_params / 1000000000.0,
                )
                return "int8"
        except Exception:
            pass
    return "none"


def _build_bnb_config(quantization: str):
    if quantization not in ("int8", "int4"):
        return None
    try:
        from transformers import BitsAndBytesConfig
        import torch
    except ImportError as exc:
        raise RuntimeError(
            "量化需要 bitsandbytes 和 transformers\n请安装：pip install bitsandbytes transformers accelerate"
        ) from exc
    if quantization == "int8":
        return BitsAndBytesConfig(load_in_8bit=True)
    return BitsAndBytesConfig(
        load_in_4bit=True,
        bnb_4bit_compute_dtype=torch.float16,
        bnb_4bit_use_double_quant=True,
        bnb_4bit_quant_type="nf4",
    )


def _log_vram_estimate(quantization: str, model_name: str) -> None:
    estimates = {
        ("7b", "none"): "~14 GB (fp16)",
        (
            "20b",
            "none",
        ): "~38 GB (fp16) — 单张 H100 可用；L40S 也可装下，但 FP8 会降级至 BF16",
        ("20b", "int8"): "~19 GB (int8) — L40S 充裕，精度损失 <1%%",
        ("40b", "none"): "~80 GB (fp16) — L40S 单卡装不下！",
        ("40b", "int8"): "~40 GB (int8) — L40S 勉强可用，建议预留 4GB overhead",
        ("40b", "int4"): "~22 GB (int4/NF4) — L40S 充裕，可与 Boltz2 并行",
    }
    if "40b" in model_name.lower():
        size_key = "40b"
    elif "20b" in model_name.lower():
        size_key = "20b"
    else:
        size_key = "7b"
    est = estimates.get((size_key, quantization), "未知")
    logger.info("[显存预估] %s + %s → %s", model_name, quantization, est)


def _resolve_evo2_model_name(model_dir: Path) -> str:
    valid_names = (
        "evo2_40b",
        "evo2_7b",
        "evo2_20b",
        "evo2_40b_base",
        "evo2_7b_base",
        "evo2_1b_base",
        "evo2_7b_262k",
        "evo2_7b_microviridae",
    )
    if model_dir.name in valid_names:
        return model_dir.name
    for checkpoint in sorted(model_dir.glob("*.pt")):
        if checkpoint.stem in valid_names:
            return checkpoint.stem
    dir_name = model_dir.name.lower()
    if "40b" in dir_name:
        return "evo2_40b"
    if "20b" in dir_name:
        return "evo2_20b"
    return "evo2_7b"


def _find_local_evo2_checkpoint(model_dir: Path, model_name: str) -> Optional[Path]:
    preferred = model_dir / f"{model_name}.pt"
    if preferred.exists():
        return preferred
    checkpoints = sorted(model_dir.glob("*.pt"))
    return checkpoints[0] if checkpoints else None


def _ensure_vortex_flash_attn_shim() -> bool:
    for site_dir in site.getsitepackages():
        base = Path(site_dir)
        for candidate in sorted(
            (base / "vortex" / "ops" / "depr_attn").glob("flash_attn_2_cuda*.so")
        ):
            target = base / candidate.name
            try:
                if not target.exists():
                    try:
                        target.symlink_to(candidate)
                    except OSError:
                        shutil.copy2(candidate, target)
                return True
            except OSError:
                continue
    return False


def _load_model(model_dir: Path, quantization: str = "auto") -> None:
    global _model, _tokenizer, _model_name, _device, _quantization
    import torch

    if torch.cuda.is_available():
        _device = "cuda"
        gpu_name = torch.cuda.get_device_name(0)
        vram_gb = torch.cuda.get_device_properties(0).total_memory / 1024**3
        logger.info("[GPU] %s  VRAM: %.1fGB", gpu_name, vram_gb)
    else:
        _device = "cpu"
        logger.warning("[GPU] CUDA 不可用，使用 CPU（推理速度将极慢）")
        if quantization in ("int8", "int4"):
            logger.warning("[量化] bitsandbytes 量化需要 CUDA，已回退为 none")
            quantization = "none"
    _quantization = _detect_quantization(model_dir, quantization)
    _model_name = _resolve_evo2_model_name(model_dir)
    _log_vram_estimate(_quantization, _model_name)
    local_checkpoint = _find_local_evo2_checkpoint(model_dir, _model_name)
    evo2_import_error: Optional[Exception] = None
    Evo2 = None
    try:
        from evo2 import Evo2 as ImportedEvo2

        Evo2 = ImportedEvo2
    except (ImportError, ModuleNotFoundError) as exc:
        evo2_import_error = exc
        if "flash_attn_2_cuda" in str(exc) and _ensure_vortex_flash_attn_shim():
            importlib.invalidate_caches()
            try:
                from evo2 import Evo2 as ImportedEvo2

                Evo2 = ImportedEvo2
                logger.info(
                    "[Evo2] 已暴露 vortex 自带 flash_attn_2_cuda 扩展，重试导入成功"
                )
                evo2_import_error = None
            except Exception as retry_exc:
                evo2_import_error = retry_exc
    try:
        if Evo2 is None and evo2_import_error is not None:
            raise evo2_import_error
        if local_checkpoint is not None:
            logger.info(
                "[Evo2] 官方包加载：model=%s, checkpoint=%s（量化=%s）",
                _model_name,
                local_checkpoint,
                _quantization,
            )
        else:
            logger.info(
                "[Evo2] 官方包加载：model=%s（量化=%s）", _model_name, _quantization
            )
        if _quantization == "none":
            if local_checkpoint is not None:
                _model = Evo2(_model_name, local_path=str(local_checkpoint))
            else:
                _model = Evo2(_model_name)
        else:
            try:
                _model = Evo2(_model_name, quantization=_quantization)
            except TypeError:
                logger.warning(
                    "[Evo2] 官方包不支持 quantization= 参数，转为 transformers 后端加载量化模型"
                )
                raise ImportError("use_transformers_for_quant")
        logger.info("[Evo2] 模型加载完成（evo2 官方包）")
        return
    except (ImportError, ModuleNotFoundError) as exc:
        if local_checkpoint is not None:
            if "libcudart.so.12" in str(exc):
                raise RuntimeError(
                    "检测到本地 Evo2 checkpoint，但当前运行时缺少 CUDA 12 动态库。请安装 CUDA 12 系列的 PyTorch/运行时，或重新执行 scripts/setup_server.sh。"
                ) from exc
            raise RuntimeError(
                f"检测到本地 Evo2 checkpoint（*.pt），但 evo2 官方加载器不可用：{exc}。请检查 flash-attn / evo2 运行时依赖。"
            ) from exc
        logger.warning(
            "[Evo2] evo2 包未安装或需要 transformers 量化后端，切换备选方案..."
        )
    try:
        from transformers import AutoModelForCausalLM, AutoTokenizer
        import torch

        model_path = (
            str(model_dir)
            if model_dir.exists()
            else os.getenv("EVO2_MODEL_REPO", "arcinstitute/evo2-7b")
        )
        logger.info(
            "[Evo2] transformers 加载：%s（量化=%s）", model_path, _quantization
        )
        _tokenizer = AutoTokenizer.from_pretrained(
            model_path, trust_remote_code=True, use_fast=False
        )
        bnb_config = _build_bnb_config(_quantization)
        if bnb_config is not None:
            _model = AutoModelForCausalLM.from_pretrained(
                model_path,
                quantization_config=bnb_config,
                device_map="auto",
                trust_remote_code=True,
            )
        else:
            _model = AutoModelForCausalLM.from_pretrained(
                model_path,
                torch_dtype=torch.float16,
                device_map="auto",
                trust_remote_code=True,
            )
        _model.eval()
        logger.info("[Evo2] 模型加载完成（transformers，量化=%s）", _quantization)
        return
    except Exception as exc:
        logger.error("[Evo2] 模型加载失败：%s", exc)
        raise RuntimeError(f"Evo2 模型加载失败，请检查安装和权重路径：{exc}") from exc


def _compute_log_likelihood_evo2_pkg(sequence: str) -> float:
    scores = _model.score_sequences([sequence], batch_size=1, reduce_method="sum")
    return float(scores[0])


def _compute_log_likelihood_transformers(sequence: str) -> float:
    import torch

    inputs = _tokenizer(sequence, return_tensors="pt").to(_device)
    input_ids = inputs["input_ids"]
    with torch.no_grad():
        outputs = _model(**inputs, labels=input_ids)
    n_tokens = input_ids.shape[1] - 1
    log_likelihood = -outputs.loss.item() * n_tokens
    return float(log_likelihood)


def compute_log_likelihood(sequence: str) -> float:
    if _model is None:
        raise RuntimeError("模型尚未加载")
    model_cls = type(_model).__name__
    if "Evo2" in model_cls or hasattr(_model, "tokenize"):
        return _compute_log_likelihood_evo2_pkg(sequence)
    else:
        return _compute_log_likelihood_transformers(sequence)


def build_app(model_dir: Path, quantization: str = "auto"):
    try:
        from fastapi import FastAPI, HTTPException
    except ImportError:
        logger.error(
            "缺少 fastapi/pydantic，请安装：pip install fastapi uvicorn[standard]"
        )
        sys.exit(1)
    app = FastAPI(title="Evo2 Likelihood Server", version="1.0")
    _startup_done = {"ok": False, "error": None}

    @app.on_event("startup")
    async def _startup():
        try:
            _load_model(model_dir, quantization=quantization)
            _startup_done["ok"] = True
            logger.info(
                "[Evo2 Server] 就绪，模型=%s，量化=%s", _model_name, _quantization
            )
        except Exception as exc:
            _startup_done["error"] = str(exc)
            logger.error("[Evo2 Server] 启动失败：%s", exc)

    @app.get("/health")
    async def health():
        if _startup_done["error"]:
            raise HTTPException(status_code=503, detail=_startup_done["error"])
        if not _startup_done["ok"]:
            raise HTTPException(status_code=503, detail="模型加载中，请稍候")
        return {"status": "ok", "model": _model_name, "quantization": _quantization}

    @app.get("/ready")
    async def ready():
        return await health()

    def _score(req: LikelihoodRequest) -> LikelihoodResponse:
        if not req.sequence or not req.sequence.strip():
            raise HTTPException(status_code=422, detail="sequence 不能为空")
        seq = req.sequence.upper().strip()
        try:
            ll = compute_log_likelihood(seq)
        except Exception as exc:
            logger.exception("打分失败: %s", exc)
            raise HTTPException(status_code=500, detail=str(exc))
        return LikelihoodResponse(log_likelihood=ll, model=_model_name)

    app.post("/likelihood", response_model=LikelihoodResponse)(_score)
    app.post("/v1/likelihood", response_model=LikelihoodResponse)(_score)
    app.post("/score", response_model=LikelihoodResponse)(_score)
    app.post("/v1/score", response_model=LikelihoodResponse)(_score)
    return app


def main() -> None:
    parser = argparse.ArgumentParser(description="Evo2 推理 HTTP 服务")
    parser.add_argument(
        "--model-dir",
        type=Path,
        default=Path(os.getenv("MODEL_DIR", "/root/gpufree-data/models/evo2_7b")),
        help="模型权重目录",
    )
    parser.add_argument(
        "--port", type=int, default=int(os.getenv("MODEL_PORT", "8000"))
    )
    parser.add_argument("--host", default="0.0.0.0")
    parser.add_argument(
        "--workers", type=int, default=1, help="Uvicorn worker 数（GPU 服务建议 1）"
    )
    parser.add_argument(
        "--quantization",
        choices=["auto", "none", "int8", "int4"],
        default=os.getenv("EVO2_QUANTIZATION", "auto"),
        help="量化模式：auto=自动检测（40B 目录自动 int8）; none=fp16 全精度; int8=LLM.int8()（40B≈40GB）; int4=NF4双量化（40B≈22GB）",
    )
    args = parser.parse_args()
    logger.info("[Evo2 Server] 模型目录：%s", args.model_dir)
    logger.info("[Evo2 Server] 量化模式：%s", args.quantization)
    logger.info("[Evo2 Server] 监听：%s:%d", args.host, args.port)
    vram_limit = os.getenv("GPU_MEMORY_LIMIT_GB", "44")
    logger.info("[Evo2 Server] 显存软限制：%s GB", vram_limit)
    try:
        import uvicorn
    except ImportError:
        logger.error("缺少 uvicorn，请安装：pip install uvicorn[standard]")
        sys.exit(1)
    app = build_app(args.model_dir, quantization=args.quantization)
    uvicorn.run(
        app, host=args.host, port=args.port, workers=args.workers, log_level="info"
    )


if __name__ == "__main__":
    main()
