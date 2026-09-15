from __future__ import annotations
import csv
import json
import os
import re
import shutil
import signal
import sqlite3
import subprocess
import sys
from datetime import datetime
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any, Dict, List, Optional
from urllib.parse import parse_qs, urlparse
import yaml

ROOT_DIR = Path(__file__).resolve().parent
CONFIG_PATH = ROOT_DIR / "config" / "config.yaml"
STATIC_DIR = ROOT_DIR / "web" / "dashboard"
RUNTIME_DIR = ROOT_DIR / ".dashboard_runtime"
STATE_PATH = RUNTIME_DIR / "pipeline_state.json"
PROJECT_PYTHON = ROOT_DIR / ".venv" / "bin" / "python"
DEFAULT_HOST = "127.0.0.1"
DEFAULT_PORT = 8787
DEFAULT_LOG_LINES = 200
MAX_LOG_LINES = 2000
MAX_QUERY_ROWS = 200
ALLOWED_STAGE_VALUES = {"all", "evo", "boltz", "analyze"}
ALLOWED_EVO2_SIZES = {"7b", "20b", "40b"}
ALLOWED_EVO2_QUANTIZATION = {"auto", "none", "int8", "int4"}


def load_yaml_config(config_path: Path = CONFIG_PATH) -> Dict[str, Any]:
    with config_path.open("r", encoding="utf-8") as fp:
        return yaml.safe_load(fp) or {}


def get_pipeline_python() -> str:
    if PROJECT_PYTHON.exists():
        return str(PROJECT_PYTHON)
    return sys.executable


def ensure_runtime_dir() -> None:
    RUNTIME_DIR.mkdir(parents=True, exist_ok=True)


def now_iso() -> str:
    return datetime.now().isoformat(timespec="seconds")


def read_json_file(path: Path) -> Optional[Dict[str, Any]]:
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None


def write_json_file(path: Path, payload: Dict[str, Any]) -> None:
    ensure_runtime_dir()
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def remove_state_file() -> None:
    if STATE_PATH.exists():
        STATE_PATH.unlink()


def process_alive(pid: int) -> bool:
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    return True


def detect_external_pipeline() -> Optional[Dict[str, Any]]:
    try:
        result = run_command(["ps", "-eo", "pid=,args="])
    except Exception:
        return None
    if result.returncode != 0:
        return None
    for line in result.stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        match = re.match("^(\\d+)\\s+(.*)$", line)
        if not match:
            continue
        pid = int(match.group(1))
        args = match.group(2)
        if "run_pipeline.py" not in args:
            continue
        if "dashboard_server.py" in args:
            continue
        if pid == os.getpid():
            continue
        stage_match = re.search("--stage\\s+(\\w+)", args)
        limit_match = re.search("--limit\\s+(\\d+)", args)
        size_match = re.search("--evo2-size\\s+(\\w+)", args)
        return {
            "pid": pid,
            "running": True,
            "alive": True,
            "stage": stage_match.group(1) if stage_match else "unknown",
            "limit": int(limit_match.group(1)) if limit_match else None,
            "shutdown_on_complete": "--shutdown-on-complete" in args,
            "evo2_size": size_match.group(1) if size_match else None,
            "command_display": args,
            "source": "external-detected",
        }
    return None


def read_pipeline_state() -> Dict[str, Any]:
    payload = read_json_file(STATE_PATH) or {}
    pid = int(payload.get("pid") or 0)
    if pid and (not process_alive(pid)):
        payload["running"] = False
        payload["ended_at"] = payload.get("ended_at") or now_iso()
        write_json_file(STATE_PATH, payload)
    if not payload.get("running"):
        external = detect_external_pipeline()
        if external:
            return external
    return payload


def list_log_files() -> List[Dict[str, Any]]:
    log_dir = ROOT_DIR / "logs"
    entries: List[Dict[str, Any]] = []
    if not log_dir.exists():
        return entries
    for path in sorted(
        log_dir.glob("*.log"), key=lambda item: item.stat().st_mtime, reverse=True
    ):
        stat = path.stat()
        entries.append(
            {
                "name": path.name,
                "path": str(path),
                "size_bytes": stat.st_size,
                "modified_at": datetime.fromtimestamp(stat.st_mtime).isoformat(
                    timespec="seconds"
                ),
            }
        )
    return entries


def clamp_int(value: Optional[str], default: int, minimum: int, maximum: int) -> int:
    try:
        parsed = int(value) if value is not None else default
    except (TypeError, ValueError):
        parsed = default
    return max(minimum, min(parsed, maximum))


def tail_lines(path: Path, line_count: int) -> List[str]:
    if not path.exists():
        return []
    text = path.read_text(encoding="utf-8", errors="replace")
    return text.splitlines()[-line_count:]


def read_metrics_events(metrics_path: Path, limit: int = 50) -> List[Dict[str, Any]]:
    if not metrics_path.exists():
        return []
    events: List[Dict[str, Any]] = []
    for line in metrics_path.read_text(encoding="utf-8", errors="replace").splitlines()[
        -limit:
    ]:
        line = line.strip()
        if not line:
            continue
        try:
            events.append(json.loads(line))
        except json.JSONDecodeError:
            continue
    return events


def run_command(
    command: List[str], timeout: int = 10
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command, capture_output=True, text=True, timeout=timeout, check=False
    )


def parse_gpu_status() -> Dict[str, Any]:
    query = "index,name,memory.total,memory.free,memory.used,utilization.gpu,temperature.gpu"
    try:
        result = run_command(
            ["nvidia-smi", f"--query-gpu={query}", "--format=csv,noheader,nounits"]
        )
    except FileNotFoundError:
        return {"available": False, "error": "nvidia-smi not found", "gpus": []}
    except subprocess.TimeoutExpired:
        return {"available": False, "error": "nvidia-smi timeout", "gpus": []}
    if result.returncode != 0:
        return {
            "available": False,
            "error": result.stderr.strip() or "nvidia-smi failed",
            "gpus": [],
        }
    gpus: List[Dict[str, Any]] = []
    for line in result.stdout.splitlines():
        parts = [part.strip() for part in line.split(",")]
        if len(parts) != 7:
            continue
        total = float(parts[2])
        free = float(parts[3])
        used = float(parts[4])
        gpus.append(
            {
                "index": parts[0],
                "name": parts[1],
                "memory_total_mb": total,
                "memory_free_mb": free,
                "memory_used_mb": used,
                "memory_used_percent": round(used / total * 100.0, 2) if total else 0.0,
                "utilization_gpu_percent": float(parts[5]),
                "temperature_c": float(parts[6]),
            }
        )
    return {"available": True, "gpus": gpus}


def get_disk_usage() -> Dict[str, Any]:
    shares = {
        "code": ROOT_DIR,
        "share_data": Path("/root/gpufree-share/amr_hunter_workspace/amr_hunter_data"),
        "fast_data": Path("/root/gpufree-data"),
    }
    payload: Dict[str, Any] = {}
    for name, path in shares.items():
        if not path.exists():
            payload[name] = {"path": str(path), "exists": False}
            continue
        usage = shutil.disk_usage(path)
        payload[name] = {
            "path": str(path),
            "exists": True,
            "total_gb": round(usage.total / 1024**3, 2),
            "used_gb": round(usage.used / 1024**3, 2),
            "free_gb": round(usage.free / 1024**3, 2),
        }
    return payload


def get_db_connection(db_path: Path) -> sqlite3.Connection:
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    return conn


def summarize_database(config: Dict[str, Any]) -> Dict[str, Any]:
    db_path = Path(config["paths"]["share_db_path"])
    if not db_path.exists():
        return {"exists": False, "path": str(db_path)}
    conn = get_db_connection(db_path)
    try:
        status_counts = [
            dict(row)
            for row in conn.execute(
                "SELECT status, COUNT(*) AS count FROM mutations GROUP BY status ORDER BY count DESC"
            ).fetchall()
        ]
        gene_count = conn.execute("SELECT COUNT(*) FROM genes").fetchone()[0]
        top_completed_genes = [
            dict(row)
            for row in conn.execute(
                "\n                SELECT g.name AS gene_name, COUNT(*) AS completed_count\n                FROM mutations m\n                JOIN genes g ON g.id = m.gene_id\n                WHERE m.status = 'COMPLETED'\n                GROUP BY g.name\n                ORDER BY completed_count DESC, g.name ASC\n                LIMIT 10\n                "
            ).fetchall()
        ]
    finally:
        conn.close()
    counts_by_status = {row["status"]: row["count"] for row in status_counts}
    summary_row = {
        "mutation_count": sum(counts_by_status.values()),
        "gene_count": gene_count,
        "boltz_core_count": counts_by_status.get("COMPLETED", 0),
        "boltz_ready_count": counts_by_status.get("BOLTZ_READY", 0),
        "completed_count": counts_by_status.get("COMPLETED", 0),
        "failed_count": counts_by_status.get("FAILED", 0),
        "last_updated_at": None,
    }
    return {
        "exists": True,
        "path": str(db_path),
        "status_counts": status_counts,
        "summary": summary_row,
        "top_ready_genes": [],
        "top_completed_genes": top_completed_genes,
    }


def query_mutations(
    config: Dict[str, Any],
    *,
    status: str = "",
    gene: str = "",
    limit: int = 100,
    offset: int = 0,
) -> Dict[str, Any]:
    db_path = Path(config["paths"]["share_db_path"])
    if not db_path.exists():
        return {"rows": [], "total": 0, "path": str(db_path), "exists": False}
    conditions: List[str] = []
    params: List[Any] = []
    if status:
        conditions.append("m.status = ?")
        params.append(status)
    if gene:
        conditions.append("g.name LIKE ?")
        params.append(f"%{gene}%")
    where_sql = f"WHERE {' AND '.join(conditions)}" if conditions else ""
    conn = get_db_connection(db_path)
    try:
        total = conn.execute(
            f"SELECT COUNT(*) FROM mutations m JOIN genes g ON g.id = m.gene_id {where_sql}",
            params,
        ).fetchone()[0]
        rows = [
            dict(row)
            for row in conn.execute(
                f"\n                SELECT m.id, g.name AS gene_name, m.aa_change, m.region_type, m.region_name,\n                       m.pos, m.ref, m.alt, m.status, m.evo_model, m.evo_delta,\n                       m.boltz_ptm, m.boltz_iptm, m.boltz_plddt, m.final_interpretation,\n                       m.updated_at, m.last_error\n                FROM mutations m\n                JOIN genes g ON g.id = m.gene_id\n                {where_sql}\n                ORDER BY m.updated_at DESC, m.id DESC\n                LIMIT ? OFFSET ?\n                ",
                [*params, limit, offset],
            ).fetchall()
        ]
    finally:
        conn.close()
    return {
        "rows": rows,
        "total": total,
        "limit": limit,
        "offset": offset,
        "exists": True,
        "path": str(db_path),
    }


def execute_readonly_query(config: Dict[str, Any], sql: str) -> Dict[str, Any]:
    normalized = sql.strip().lower()
    if not normalized:
        raise ValueError("SQL cannot be empty")
    if ";" in normalized.rstrip(";"):
        raise ValueError("Only a single statement is allowed")
    if not normalized.startswith(("select", "with", "pragma", "explain")):
        raise ValueError(
            "Only read-only SELECT/WITH/PRAGMA/EXPLAIN queries are allowed"
        )
    db_path = Path(config["paths"]["share_db_path"])
    conn = get_db_connection(db_path)
    try:
        cur = conn.execute(sql)
        columns = [item[0] for item in cur.description] if cur.description else []
        rows = [dict(row) for row in cur.fetchmany(MAX_QUERY_ROWS)]
    finally:
        conn.close()
    return {
        "columns": columns,
        "rows": rows,
        "row_count": len(rows),
        "truncated": len(rows) >= MAX_QUERY_ROWS,
    }


def summarize_results(config: Dict[str, Any]) -> Dict[str, Any]:
    pipeline_cfg = config.get("pipeline") or {}
    export_csv_path = Path(str(pipeline_cfg.get("export_csv_path") or ""))
    metrics_path = Path(
        str(
            pipeline_cfg.get("metrics_path")
            or ROOT_DIR / "logs" / "batch_metrics.jsonl"
        )
    )
    artifacts_dir = (
        export_csv_path.parent / "boltz_artifacts"
        if export_csv_path.parent
        else ROOT_DIR / "logs"
    )
    summary: Dict[str, Any] = {
        "export_csv": {
            "path": str(export_csv_path),
            "exists": export_csv_path.exists(),
        },
        "metrics_path": {
            "path": str(metrics_path),
            "exists": metrics_path.exists(),
            "recent_events": read_metrics_events(metrics_path, 40),
        },
        "artifacts": {
            "path": str(artifacts_dir),
            "exists": artifacts_dir.exists(),
            "count": 0,
        },
        "csv_preview": [],
    }
    if export_csv_path.exists():
        stat = export_csv_path.stat()
        summary["export_csv"].update(
            {
                "size_bytes": stat.st_size,
                "modified_at": datetime.fromtimestamp(stat.st_mtime).isoformat(
                    timespec="seconds"
                ),
            }
        )
        with export_csv_path.open("r", encoding="utf-8", errors="replace") as fp:
            reader = csv.DictReader(fp)
            summary["csv_preview"] = [row for _, row in zip(range(10), reader)]
    if artifacts_dir.exists():
        summary["artifacts"]["count"] = sum((1 for _ in artifacts_dir.iterdir()))
    return summary


def start_pipeline(payload: Dict[str, Any]) -> Dict[str, Any]:
    current = read_pipeline_state()
    if current.get("running") and process_alive(int(current.get("pid") or 0)):
        raise RuntimeError("A pipeline task is already running")
    stage = str(payload.get("stage") or "boltz")
    if stage not in ALLOWED_STAGE_VALUES:
        raise ValueError(f"Unsupported stage: {stage}")
    evo2_size = str(payload.get("evo2_size") or "7b")
    if evo2_size not in ALLOWED_EVO2_SIZES:
        raise ValueError(f"Unsupported evo2_size: {evo2_size}")
    evo2_quantization = str(payload.get("evo2_quantization") or "auto")
    if evo2_quantization not in ALLOWED_EVO2_QUANTIZATION:
        raise ValueError(f"Unsupported evo2_quantization: {evo2_quantization}")
    limit = payload.get("limit")
    if limit in ("", None):
        limit = None
    else:
        limit = int(limit)
        if limit <= 0:
            raise ValueError("limit must be > 0")
    shutdown_on_complete = bool(payload.get("shutdown_on_complete", False))
    log_name = str(
        payload.get("log_name")
        or f"dashboard_{stage}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    )
    if "/" in log_name or ".." in log_name:
        raise ValueError("Invalid log_name")
    log_path = ROOT_DIR / "logs" / log_name
    log_path.parent.mkdir(parents=True, exist_ok=True)
    cmd: List[str] = [
        get_pipeline_python(),
        str(ROOT_DIR / "run_pipeline.py"),
        "--stage",
        stage,
        "--evo2-size",
        evo2_size,
    ]
    if evo2_quantization:
        cmd.extend(["--evo2-quantization", evo2_quantization])
    if limit is not None:
        cmd.extend(["--limit", str(limit)])
    if shutdown_on_complete:
        cmd.append("--shutdown-on-complete")
    stdout_handle = log_path.open("a", encoding="utf-8")
    process = subprocess.Popen(
        cmd,
        cwd=ROOT_DIR,
        stdout=stdout_handle,
        stderr=subprocess.STDOUT,
        start_new_session=True,
    )
    stdout_handle.close()
    state = {
        "pid": process.pid,
        "running": True,
        "command": cmd,
        "command_display": " ".join(cmd),
        "stage": stage,
        "limit": limit,
        "shutdown_on_complete": shutdown_on_complete,
        "evo2_size": evo2_size,
        "evo2_quantization": evo2_quantization,
        "log_name": log_name,
        "log_path": str(log_path),
        "started_at": now_iso(),
    }
    write_json_file(STATE_PATH, state)
    return state


def stop_pipeline() -> Dict[str, Any]:
    state = read_pipeline_state()
    pid = int(state.get("pid") or 0)
    if not pid or not process_alive(pid):
        state["running"] = False
        write_json_file(STATE_PATH, state)
        raise RuntimeError("No running pipeline process found")
    os.killpg(pid, signal.SIGTERM)
    state["running"] = False
    state["stop_requested_at"] = now_iso()
    write_json_file(STATE_PATH, state)
    return state


def build_overview() -> Dict[str, Any]:
    config = load_yaml_config()
    pipeline_state = read_pipeline_state()
    pipeline_state["alive"] = process_alive(int(pipeline_state.get("pid") or 0))
    return {
        "server_time": now_iso(),
        "pipeline": pipeline_state,
        "gpu": parse_gpu_status(),
        "disk": get_disk_usage(),
        "database": summarize_database(config),
        "results": summarize_results(config),
        "logs": list_log_files(),
    }


class DashboardHandler(BaseHTTPRequestHandler):
    server_version = "AMRHunterDashboard/0.1"

    def log_message(self, format: str, *args: Any) -> None:
        return

    def _json_response(self, payload: Any, status: int = HTTPStatus.OK) -> None:
        body = json.dumps(payload, ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _text_response(
        self,
        payload: str,
        *,
        content_type: str = "text/html; charset=utf-8",
        status: int = HTTPStatus.OK,
    ) -> None:
        body = payload.encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _read_json_body(self) -> Dict[str, Any]:
        content_length = int(self.headers.get("Content-Length") or 0)
        if content_length <= 0:
            return {}
        raw = self.rfile.read(content_length)
        if not raw:
            return {}
        return json.loads(raw.decode("utf-8"))

    def _serve_static(self, relative_path: str) -> None:
        path = (STATIC_DIR / relative_path).resolve()
        if (
            not str(path).startswith(str(STATIC_DIR.resolve()))
            or not path.exists()
            or (not path.is_file())
        ):
            self._json_response({"error": "Not found"}, status=HTTPStatus.NOT_FOUND)
            return
        content_type = "text/plain; charset=utf-8"
        if path.suffix == ".html":
            content_type = "text/html; charset=utf-8"
        elif path.suffix == ".css":
            content_type = "text/css; charset=utf-8"
        elif path.suffix == ".js":
            content_type = "application/javascript; charset=utf-8"
        self._text_response(path.read_text(encoding="utf-8"), content_type=content_type)

    def do_GET(self) -> None:
        parsed = urlparse(self.path)
        query = parse_qs(parsed.query)
        try:
            if parsed.path == "/":
                self._serve_static("index.html")
                return
            if parsed.path.startswith("/static/"):
                self._serve_static(parsed.path.removeprefix("/static/"))
                return
            if parsed.path == "/api/overview":
                self._json_response(build_overview())
                return
            if parsed.path == "/api/gpu":
                self._json_response(parse_gpu_status())
                return
            if parsed.path == "/api/logs":
                log_name = str(query.get("name", [""])[0])
                if not log_name or "/" in log_name or ".." in log_name:
                    raise ValueError("A valid log file name is required")
                line_count = clamp_int(
                    query.get("lines", [None])[0], DEFAULT_LOG_LINES, 10, MAX_LOG_LINES
                )
                log_path = ROOT_DIR / "logs" / log_name
                self._json_response(
                    {
                        "name": log_name,
                        "path": str(log_path),
                        "exists": log_path.exists(),
                        "lines": tail_lines(log_path, line_count),
                    }
                )
                return
            if parsed.path == "/api/database/mutations":
                config = load_yaml_config()
                status = str(query.get("status", [""])[0])
                gene = str(query.get("gene", [""])[0])
                limit = clamp_int(query.get("limit", [None])[0], 100, 1, 500)
                offset = clamp_int(query.get("offset", [None])[0], 0, 0, 100000)
                self._json_response(
                    query_mutations(
                        config, status=status, gene=gene, limit=limit, offset=offset
                    )
                )
                return
            if parsed.path == "/api/results":
                self._json_response(summarize_results(load_yaml_config()))
                return
            self._json_response({"error": "Not found"}, status=HTTPStatus.NOT_FOUND)
        except Exception as exc:
            self._json_response({"error": str(exc)}, status=HTTPStatus.BAD_REQUEST)

    def do_POST(self) -> None:
        parsed = urlparse(self.path)
        try:
            payload = self._read_json_body()
            if parsed.path == "/api/pipeline/start":
                self._json_response(start_pipeline(payload), status=HTTPStatus.CREATED)
                return
            if parsed.path == "/api/pipeline/stop":
                self._json_response(stop_pipeline())
                return
            if parsed.path == "/api/database/query":
                sql = str(payload.get("sql") or "")
                self._json_response(execute_readonly_query(load_yaml_config(), sql))
                return
            self._json_response({"error": "Not found"}, status=HTTPStatus.NOT_FOUND)
        except RuntimeError as exc:
            self._json_response({"error": str(exc)}, status=HTTPStatus.CONFLICT)
        except Exception as exc:
            self._json_response({"error": str(exc)}, status=HTTPStatus.BAD_REQUEST)


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser(description="AMR-Hunter 项目管理控制台")
    parser.add_argument("--host", default=DEFAULT_HOST, help="监听地址，默认 127.0.0.1")
    parser.add_argument(
        "--port", type=int, default=DEFAULT_PORT, help="监听端口，默认 8787"
    )
    args = parser.parse_args()
    ensure_runtime_dir()
    server = ThreadingHTTPServer((args.host, args.port), DashboardHandler)
    print(f"AMR-Hunter dashboard listening on http://{args.host}:{args.port}")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
