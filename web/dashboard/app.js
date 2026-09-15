const state = {
  overview: null,
};

function $(id) {
  return document.getElementById(id);
}

function formatNumber(value) {
  if (value === null || value === undefined || value === "") {
    return "--";
  }
  if (typeof value === "number") {
    return Number.isInteger(value) ? value.toString() : value.toFixed(3);
  }
  return String(value);
}

function formatJson(value) {
  return JSON.stringify(value, null, 2);
}

async function apiFetch(url, options = {}) {
  const response = await fetch(url, {
    headers: { "Content-Type": "application/json", ...(options.headers || {}) },
    ...options,
  });
  const payload = await response.json();
  if (!response.ok) {
    throw new Error(payload.error || `Request failed: ${response.status}`);
  }
  return payload;
}

function renderSimpleCards(container, items, emptyText) {
  if (!items || items.length === 0) {
    container.innerHTML = `<div class="empty-state">${emptyText}</div>`;
    return;
  }
  container.innerHTML = items
    .map(
      (item) => `
        <div class="summary-card">
          <div class="label">${item.label}</div>
          <div class="value">${item.value}</div>
          ${item.meta ? `<div class="meta">${item.meta}</div>` : ""}
        </div>
      `,
    )
    .join("");
}

function renderPipelineSummary(pipeline) {
  const container = $("pipeline-summary");
  const badge = $("pipeline-badge");
  const alive = Boolean(pipeline && pipeline.alive);
  badge.textContent = alive ? "运行中" : pipeline?.pid ? "已退出" : "未运行";
  badge.className = `status-badge ${alive ? "running" : pipeline?.pid ? "stopped" : "idle"}`;

  const items = [
    { label: "PID", value: pipeline?.pid || "--" },
    { label: "Stage", value: pipeline?.stage || "--" },
    { label: "Limit", value: pipeline?.limit || "跑空队列" },
    { label: "自动关机", value: pipeline?.shutdown_on_complete ? "开启" : "关闭" },
    { label: "开始时间", value: pipeline?.started_at || "--" },
    { label: "日志", value: pipeline?.log_name || "--" },
  ];
  renderSimpleCards(container, items, "暂无运行记录");
}

function renderGpu(gpu) {
  const container = $("gpu-cards");
  if (!gpu.available) {
    container.innerHTML = `<div class="empty-state">GPU 状态不可用：${gpu.error || "unknown"}</div>`;
    return;
  }
  container.innerHTML = gpu.gpus
    .map(
      (item) => `
        <div class="gpu-card">
          <div class="label">GPU ${item.index}</div>
          <div class="value">${item.name}</div>
          <div class="meta">显存 ${item.memory_used_mb}/${item.memory_total_mb} MB (${item.memory_used_percent}%)</div>
          <div class="meta">利用率 ${item.utilization_gpu_percent}% / 温度 ${item.temperature_c} C</div>
        </div>
      `,
    )
    .join("");
}

function renderDisk(disk) {
  const container = $("disk-cards");
  const items = Object.entries(disk || {}).map(([name, info]) => ({
    label: name,
    value: info.exists ? `${info.free_gb} GB free` : "missing",
    meta: info.path,
  }));
  renderSimpleCards(container, items, "无磁盘信息");
}

function renderStatusCounts(counts) {
  const container = $("status-counts");
  if (!counts || counts.length === 0) {
    container.innerHTML = '<div class="empty-state">暂无状态数据</div>';
    return;
  }
  container.innerHTML = counts
    .map(
      (item) => `
        <div class="list-item">
          <div>
            <strong>${item.status}</strong>
          </div>
          <div>${item.count}</div>
        </div>
      `,
    )
    .join("");
}

function renderNamedCounts(containerId, rows, valueField) {
  const container = $(containerId);
  if (!rows || rows.length === 0) {
    container.innerHTML = '<div class="empty-state">暂无数据</div>';
    return;
  }
  container.innerHTML = rows
    .map(
      (row) => `
        <div class="list-item">
          <div>
            <strong>${row.gene_name}</strong>
          </div>
          <div>${row[valueField]}</div>
        </div>
      `,
    )
    .join("");
}

function renderKpis(database, results) {
  const summary = database.summary || {};
  const kpis = [
    { label: "总突变数", value: formatNumber(summary.mutation_count) },
    { label: "BOLTZ_READY", value: formatNumber(summary.boltz_ready_count) },
    { label: "已完成", value: formatNumber(summary.completed_count) },
    { label: "Boltz 核心结果", value: formatNumber(summary.boltz_core_count) },
    { label: "失败", value: formatNumber(summary.failed_count) },
    { label: "Boltz 产物目录", value: formatNumber(results.artifacts?.count) },
  ];
  renderSimpleCards($("result-kpis"), kpis, "暂无 KPI");
}

function renderResultMeta(results) {
  const exportCsv = results.export_csv || {};
  const meta = {
    export_csv: exportCsv,
    artifacts: results.artifacts,
    metrics_path: results.metrics_path?.path,
  };
  $("result-file-meta").textContent = formatJson(meta);
  $("metrics-events").textContent = formatJson(results.metrics_path?.recent_events || []);
}

function populateLogSelect(logs) {
  const select = $("log-select");
  const current = select.value;
  select.innerHTML = logs
    .map((item) => `<option value="${item.name}">${item.name}</option>`)
    .join("");
  if (logs.length === 0) {
    select.innerHTML = '<option value="">无日志文件</option>';
    return;
  }
  const preferred = logs.find((item) => item.name.includes("dashboard") || item.name.includes("run_pipeline"));
  select.value = logs.some((item) => item.name === current) ? current : preferred?.name || logs[0].name;
}

async function loadLog() {
  const name = $("log-select").value;
  if (!name) {
    $("log-viewer").textContent = "暂无日志可显示";
    return;
  }
  const lines = $("log-lines").value || 200;
  const payload = await apiFetch(`/api/logs?name=${encodeURIComponent(name)}&lines=${encodeURIComponent(lines)}`);
  $("log-viewer").textContent = (payload.lines || []).join("\n") || "日志为空";
}

function renderMutations(rows) {
  const tbody = $("mutation-table");
  if (!rows || rows.length === 0) {
    tbody.innerHTML = '<tr><td colspan="8">没有匹配记录</td></tr>';
    return;
  }
  tbody.innerHTML = rows
    .map(
      (row) => `
        <tr>
          <td>${row.id}</td>
          <td>${row.gene_name}</td>
          <td>${row.aa_change || `${row.ref}${row.pos}${row.alt}`}</td>
          <td>${row.status}</td>
          <td>${formatNumber(row.evo_delta)}</td>
          <td>${formatNumber(row.boltz_ptm)}</td>
          <td>${formatNumber(row.boltz_plddt)}</td>
          <td>${row.updated_at || "--"}</td>
        </tr>
      `,
    )
    .join("");
}

async function loadMutations() {
  const status = $("mutation-status").value.trim();
  const gene = $("mutation-gene").value.trim();
  const limit = $("mutation-limit").value || 100;
  const url = `/api/database/mutations?status=${encodeURIComponent(status)}&gene=${encodeURIComponent(gene)}&limit=${encodeURIComponent(limit)}`;
  const payload = await apiFetch(url);
  renderMutations(payload.rows || []);
}

async function runSql() {
  const sql = $("sql-text").value;
  const payload = await apiFetch("/api/database/query", {
    method: "POST",
    body: JSON.stringify({ sql }),
  });
  $("sql-result").textContent = formatJson(payload);
}

async function startPipeline(event) {
  event.preventDefault();
  const form = new FormData($("start-form"));
  const payload = {
    stage: form.get("stage"),
    limit: form.get("limit") || null,
    evo2_size: form.get("evo2_size"),
    evo2_quantization: form.get("evo2_quantization"),
    log_name: form.get("log_name") || null,
    shutdown_on_complete: form.get("shutdown_on_complete") === "on",
  };
  const result = await apiFetch("/api/pipeline/start", {
    method: "POST",
    body: JSON.stringify(payload),
  });
  $("control-message").textContent = `已启动：PID ${result.pid}，日志 ${result.log_name}`;
  await refreshOverview();
}

async function stopPipeline() {
  const result = await apiFetch("/api/pipeline/stop", { method: "POST", body: JSON.stringify({}) });
  $("control-message").textContent = `已发送停止信号：PID ${result.pid}`;
  await refreshOverview();
}

async function refreshOverview() {
  const payload = await apiFetch("/api/overview");
  state.overview = payload;
  $("server-time").textContent = payload.server_time;
  renderPipelineSummary(payload.pipeline || {});
  renderGpu(payload.gpu || {});
  renderDisk(payload.disk || {});
  renderKpis(payload.database || {}, payload.results || {});
  renderStatusCounts(payload.database?.status_counts || []);
  renderNamedCounts("top-ready-genes", payload.database?.top_completed_genes || [], "completed_count");
  renderResultMeta(payload.results || {});
  populateLogSelect(payload.logs || []);
  await loadLog();
  await loadMutations();
}

function wireEvents() {
  $("refresh-overview").addEventListener("click", refreshOverview);
  $("start-form").addEventListener("submit", async (event) => {
    try {
      await startPipeline(event);
    } catch (error) {
      $("control-message").textContent = error.message;
    }
  });
  $("stop-pipeline").addEventListener("click", async () => {
    try {
      await stopPipeline();
    } catch (error) {
      $("control-message").textContent = error.message;
    }
  });
  $("load-log").addEventListener("click", async () => {
    try {
      await loadLog();
    } catch (error) {
      $("log-viewer").textContent = error.message;
    }
  });
  $("load-mutations").addEventListener("click", async () => {
    try {
      await loadMutations();
    } catch (error) {
      $("mutation-table").innerHTML = `<tr><td colspan="8">${error.message}</td></tr>`;
    }
  });
  $("run-sql").addEventListener("click", async () => {
    try {
      await runSql();
    } catch (error) {
      $("sql-result").textContent = error.message;
    }
  });
}

wireEvents();
refreshOverview().catch((error) => {
  $("control-message").textContent = error.message;
});
setInterval(() => {
  refreshOverview().catch(() => {});
}, 30000);