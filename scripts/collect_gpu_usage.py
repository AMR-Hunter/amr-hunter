from __future__ import annotations
import argparse
import csv
import signal
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

RUNNING = True


def handle_signal(signum, frame):
    global RUNNING
    RUNNING = False


def run_query(args: list[str]) -> list[list[str]]:
    result = subprocess.run(args, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(result.stderr.strip() or f"command failed: {' '.join(args)}")
    rows: list[list[str]] = []
    for line in result.stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        rows.append([part.strip() for part in line.split(",")])
    return rows


def append_rows(path: Path, header: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    file_exists = path.exists() and path.stat().st_size > 0
    with path.open("a", encoding="utf-8", newline="") as fp:
        writer = csv.writer(fp)
        if not file_exists:
            writer.writerow(header)
        writer.writerows(rows)
        fp.flush()


def sample_gpu() -> list[list[str]]:
    query = "index,name,memory.used,memory.total,utilization.gpu,temperature.gpu"
    rows = run_query(
        ["nvidia-smi", f"--query-gpu={query}", "--format=csv,noheader,nounits"]
    )
    timestamp = datetime.now().isoformat(timespec="seconds")
    return [[timestamp, *row] for row in rows]


def sample_processes() -> list[list[str]]:
    query = "pid,process_name,used_memory"
    rows = run_query(
        ["nvidia-smi", f"--query-compute-apps={query}", "--format=csv,noheader,nounits"]
    )
    timestamp = datetime.now().isoformat(timespec="seconds")
    return [[timestamp, *row] for row in rows]


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Collect GPU and process memory usage into CSV files"
    )
    parser.add_argument(
        "--interval", type=int, default=30, help="Sampling interval in seconds"
    )
    parser.add_argument(
        "--gpu-output", required=True, help="CSV path for overall GPU usage"
    )
    parser.add_argument(
        "--process-output", required=True, help="CSV path for per-process GPU usage"
    )
    args = parser.parse_args()
    signal.signal(signal.SIGTERM, handle_signal)
    signal.signal(signal.SIGINT, handle_signal)
    gpu_path = Path(args.gpu_output)
    proc_path = Path(args.process_output)
    gpu_header = [
        "timestamp",
        "gpu_index",
        "gpu_name",
        "memory_used_mb",
        "memory_total_mb",
        "utilization_gpu_percent",
        "temperature_c",
    ]
    proc_header = ["timestamp", "pid", "process_name", "used_memory_mb"]
    while RUNNING:
        try:
            append_rows(gpu_path, gpu_header, sample_gpu())
            append_rows(proc_path, proc_header, sample_processes())
        except Exception as exc:
            print(f"collector error: {exc}", file=sys.stderr, flush=True)
        time.sleep(max(1, args.interval))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
