# 多阶段构建：编译依赖 + 最终运行镜像
FROM nvidia/cuda:12.1.0-base-ubuntu22.04 as builder

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3.11 \
    python3.11-venv \
    python3.11-dev \
    python3-pip \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# 创建虚拟环境（减少最终镜像大小）
RUN python3.11 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# 安装依赖（在构建阶段，减少最终镜像大小）
COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv
COPY pyproject.toml /tmp/pyproject.toml
RUN VIRTUAL_ENV=/opt/venv uv pip install -r /tmp/pyproject.toml

# ============================================================
# 最终运行镜像
FROM nvidia/cuda:12.1.0-base-ubuntu22.04

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3.11 \
    sqlite3 \
    && rm -rf /var/lib/apt/lists/*

# 复制虚拟环境
COPY --from=builder /opt/venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# 创建工作目录
WORKDIR /app

# 复制代码到容器
COPY config/ /app/config/
COPY core/ /app/core/
COPY scripts/ /app/scripts/
COPY main.py /app/main.py
COPY auto_run.py /app/auto_run.py

# 创建数据和日志目录（容器内挂载点）
RUN mkdir -p /app/data/db /app/data/results /app/data/reference /app/data/structures /app/data/ligands \
    && mkdir -p /app/logs \
    && mkdir -p /app/share_data

# 健康检查
HEALTHCHECK --interval=10s --timeout=3s --start-period=5s --retries=3 \
    CMD python3 -c "import subprocess; subprocess.run(['nvidia-smi'], check=True)" || exit 1

# 入口点
ENTRYPOINT ["python", "auto_run.py"]
CMD ["--stage", "all"]
