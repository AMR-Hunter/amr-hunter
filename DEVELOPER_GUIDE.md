# AMR-Hunter 开发者源码阅读指南

> 这份文档的目标不是“介绍功能”，而是让你**真正掌控代码**：知道每个模块做什么、数据如何流动、哪里该改、改完会影响哪里。

---

## 1. 先建立全局心智模型

把 AMR-Hunter 想成一条三段式流水线：

1) **造数据**（初始化）
- 从参考基因组提取目标基因
- 生成每个基因的单点突变
- 写入 SQLite（`PENDING`）

2) **筛数据**（Evo2）
- 对每个突变算 `evo_delta = mt_score - wt_score`
- 按 `logic_type` 做阈值分流
- 决定谁能进结构阶段（`BOLTZ_READY`）

3) **精查 + 判定**（Boltz-2 + Analyze）
- 根据 `scenario` 构造不同 payload
- 回填结构指标
- 最后做通路上位性修正 + 多基因协同修正并导出

---

## 2. 目录结构图（骨架图）

```text
amr_hunter/
├─ config/
│  └─ config.yaml
│     └─ 项目总开关：基因列表、logic_type、scenario、阈值、通路规则
│
├─ core/
│  ├─ database.py
│  │  └─ 数据库管家：建表、迁移、索引、状态约束、本地与共享库同步
│  ├─ env_manager.py
│  │  └─ 容器环境管理：路径抽象、自动探测共享路径、跨容器同步、自愈降级
│  ├─ generator.py
│  │  └─ 突变引擎：把一条 CDS 生成为全量单碱基替换突变
│  ├─ clients.py
│  │  └─ 模型客户端：Evo2/Boltz-2 HTTP 调用、重试、健康检查、Boltz payload 组装
│  └─ analyzer.py
│     └─ 业务判定器：Evo 阈值分流 + 上位性网络修正 + 最终解释
│
├─ scripts/
│  └─ init_project.py
│     └─ 初始化入口：读配置、入库基因、生成并批量写入突变
│
├─ main.py
│  └─ 主编排器：--stage evo/boltz/analyze/all 三阶段执行
│
├─ docker-compose.yml
│  └─ 模型服务编排：evo2 与 boltz2 容器（建议分阶段启动）
│
├─ pyproject.toml
│  └─ Python 依赖与版本约束（uv 管理）
│
├─ uv.lock
│  └─ uv 依赖锁文件（固定版本，保证可复现构建）
│
└─ README.md
   └─ 对外使用文档（给使用者）

├─ run_pipeline.py
│  └─ 高级编排器：管理 VRAM、处理数据库同步 (Local <=> Share)、支持 --limit、--evo2-size 参数。
│
├─ services/
│  ├─ evo2_server.py
│  │  └─ 集成了 Transformer Engine 和 FP8/BF16 动态切换的 Evo2 模型服务器。
│  └─ boltz2_server.py
│     └─ 基于 Boltz-2 的结构预测服务，集成 GPU 互斥锁。
```

---

## 3. 关键更新：模型区分与结构缓存

### 3.1 跨规格模型支持 (Evo2 7B/20B)
为了在同一份数据库中对比不同规格模型的效果，数据库架构 (v23) 引入了 `evo_model` 字段。
- 在 `mutations` 表中，由 `(mutation_id, evo_model)` 构成唯一索引。
- 支持 `7b` 和 `20b` 或其他规格的结果并存，互不覆盖。

### 3.2 结构预测缓存 (Boltz Cache)
Boltz-2 模拟非常耗时且消耗 GPU。
- **机制**：通过对 (序列 + 场景 + 配体) 生成 MD5 哈希作为 `cache_key`。
- **作用**：当不同 Evo 规格下产生相同的突变序列时（或不同基因产生相同片段时），直接从 `boltz_cache` 查表回填，无需重复运行 GPU 模拟。

---

## 4. 数据流向剖析：一个突变的一生（以 atpE 为例）

下面用“atpE 基因第 100 位 A→G”举例，完整走一遍。

---

### 阶段 A：它如何被创造出来（generator.py）

入口在 `scripts/init_project.py`：
1) 读 `config/config.yaml` 的 `target_genes`，看到 `atpE`
2) 从 `reference_genomes.path` 指向的 `.gbk` 中提取 atpE 的 CDS
3) 调用 `MutationEngine.generate_saturation_mutagenesis(...)`

`generate_saturation_mutagenesis` 的关键逻辑：
- 遍历序列每个位点
- 对每个位点尝试替换为 `A/C/G/T` 中除原碱基外的 3 种
- 计算氨基酸变化（如 `A34V`）和突变类型（同义/非同义）
- 生成记录，状态初始 `PENDING`

你可把它理解为：
- 输入一条基因序列
- 输出 N 条“候选突变任务单”

---

### 阶段 B：它如何进入数据库（database.py + init_project.py）

初始化脚本随后调用：
- `upsert_species`
- `upsert_gene`
- `MutationEngine.bulk_insert_mutations`

关键落库表：
- `genes`：存基因元信息（`logic_type`, `scenario`, `pdb_path`, `ligand_sdf_path`, `dna_sequence`）
- `mutations`：存突变任务（`pos/ref/alt/aa_change/status=PENDING`）

`database.py` 负责：
- 自动建表
- 自动迁移（旧字段兼容）
- 状态值约束（防止写入非法状态）

---

### 阶段 C：它如何被送到 Evo2（main.py + clients.py）

执行命令：
```bash
python main.py --stage evo
```

`main.py/run_evo_stage` 主要步骤：
1) 查询 `mutations.status='PENDING'`
2) 按基因分组，先计算每个基因 WT score（缓存）
3) 为每个突变构造 mutant sequence
4) 并发调用 `Evo2Client.get_likelihood(mt_seq)`
5) 回填：
   - `evo_mt_score`
   - `evo_delta = mt_score - wt_score`

之后调用 `apply_differential_filtering` 分流：
- `POSITIVE`：
  - `delta < -10.0` → `LETHAL_SKIP`
  - 否则 → `BOLTZ_READY`
- `NEGATIVE`：
  - `delta < -8.0` → `LOSS_OF_FUNCTION`（短路，直接判定耐药，不进 Boltz）
  - 否则 → `BOLTZ_READY`
- `STRUCTURAL`：
  - `delta < -10.0` → `LETHAL_SKIP`
  - `-10.0 <= delta < -8.0` → `LOSS_OF_FUNCTION`
  - 否则 → `BOLTZ_READY`

---

### 阶段 D：如何按 scenario 组装 JSON 发给 Boltz-2（clients.py）

执行命令：
```bash
python main.py --stage boltz
```

`main.py/run_boltz_stage` 只取 `BOLTZ_READY` 且 `region_type='CDS'` 的突变。

也就是说：
- 编码区突变：会进入 Boltz 结构模拟
- 非编码区（`REGULATORY`）突变：不会进入 Boltz，只保留 Evo 与分析阶段结果

对于 atpE，`config.yaml` 中 `scenario: PROTEIN_LIGAND`，因此 `Boltz2Client` 使用 `_build_ligand_payload(...)`：

```json
{
  "scenario": "protein_ligand",
  "protein": {
    "sequence": "...突变后蛋白序列...",
    "template_pdb": "/root/amr_hunter/data/structures/atpE_wt.pdb"
  },
  "ligand": {
    "sdf_path": "/root/amr_hunter/data/ligands/bedaquiline.sdf"
  }
}
```

其他场景：
- `PROTEIN_DNA` / `DIMER_DNA` → `_build_dna_payload(...)`
- `FOLDING_ONLY` / `PROTEIN_ONLY` → `_build_protein_payload(...)`

Boltz 返回后，写入：
- `boltz_iptm`
- `boltz_plddt`
- `boltz_complex_energy`
- `boltz_binding_affinity`
- `boltz_affinity_pred_value`
- `boltz_pair_energy`
- `boltz_stability_score`
- `boltz_clash`

状态更新为 `COMPLETED`。

注意：
- `boltz_binding_affinity`、`boltz_complex_energy`、`boltz_affinity_pred_value`、`boltz_pair_energy` 当前不作为核心判读字段，只保留为辅助原始观测。
- 分析导出时，这三列会以 `aux_*` 前缀输出，避免后续人工或脚本误把它们当成稳定主指标。

字段分组约定：
- 核心判读字段：`evo_delta`、`boltz_ptm`、`boltz_iptm`、`boltz_plddt`
- 专用支持字段：`boltz_stability_score`、`boltz_clash`
- 辅助留档字段：`boltz_binding_affinity`、`boltz_affinity_pred_value`、`boltz_pair_energy`、`boltz_complex_energy`、`boltz_artifact_dir`

分析阶段默认把“核心结构证据”收敛到 `boltz_ptm`、`boltz_iptm`、`boltz_plddt`；`boltz_stability_score` 只在出现明确能量型输出时，作为特定通路逻辑的补充信号。

---

### 阶段 E：最后如何判定“耐药/敏感”（analyzer.py）

执行命令：
```bash
python main.py --stage analyze
```

处理步骤：
1) 拉取 `COMPLETED / LOSS_OF_FUNCTION / LETHAL_SKIP`
2) `PathwayAnalyzer.apply_epistasis_corrections(...)` 进行通路修正
3) `PathwayAnalyzer.apply_multi_gene_synergy(...)` 进行多基因协同打分/标注
4) 特殊网络：`Rv0678 + (mmpS5/mmpL5)`
   - 如果 `Rv0678` 表现耐药（`LOSS_OF_FUNCTION`）
   - 且 mmpS5 或 mmpL5 崩溃（`LOSS_OF_FUNCTION` / `LETHAL_SKIP`）
   - 则把 Rv0678 结果翻转为敏感语义
5) 导出 CSV 前，`RE_SENSITIZED` 映射为 `SENSITIVE`

所以，一个 atpE 突变最终状态通常来自：
- Evo 阈值 + 核心结构证据（`ptm/iptm/plddt`）+ 通路修正规则（若涉及通路）
- 若存在明确能量型 `boltz_stability_score`，它只作为特定通路逻辑的补充信号，不替代核心结构证据

---

## 4. 核心数据结构：mutations 表字段详解

`mutations` 表（SQLite）关键字段如下：

- `id`：突变记录主键
- `gene_id`：关联 `genes.id`
- `region_type`：区域类型（当前约束为 `CDS` 或 `REGULATORY`）
- `region_name`：区域名称（如 `RV0678_OPERATOR_CORE`）
- `pos`：核苷酸位点（1-based）
- `ref` / `alt`：参考碱基与替换碱基
- `source_sequence`：突变所在的原始序列（用于非编码区 Evo 评分）
- `aa_change`：氨基酸变化，如 `A34V`（非编码区为空）
- `evo_mt_score`：Evo2 对突变序列评分
- `evo_delta`：`mt - wt` 差值（分流核心）
- `boltz_binding_affinity`：辅助原始亲和力指标，仅供回溯，不参与核心判读
- `boltz_affinity_pred_value`：辅助原始亲和力数值，仅供回溯，不参与核心判读
- `boltz_pair_energy`：辅助原始配对能量，仅供回溯，不参与核心判读
- `boltz_stability_score`：结构稳定性指标
- `boltz_clash`：空间冲突指标
- `boltz_iptm`：界面质量评分
- `boltz_plddt`：残基置信度评分
- `boltz_complex_energy`：辅助原始复合体能量，仅供回溯，不参与核心判读
- `synergy_score`：多基因协同得分（由协同引擎计算）
- `synergy_tags`：多基因协同标签（`|` 分隔字符串）
- `final_interpretation`：最终解释文本（可追溯）
- `target_drug_id`：关联药物（可空）
- `retry_count`：失败重试计数
- `last_error`：最后一次错误信息
- `status`：任务状态流转（最重要）
- `created_at` / `updated_at`：时间戳

状态值集合：
- `PENDING`
- `EVO_DONE`
- `BOLTZ_READY`
- `LETHAL_SKIP`
- `LOSS_OF_FUNCTION`
- `COMPLETED`
- `RE_SENSITIZED`
- `FAILED`

---

## 5. 设计模式解释（为什么这样拆）

### 5.1 为什么拆 `--stage evo` 和 `--stage boltz`

核心原因：**显存管理**。

- Evo2 和 Boltz-2 都吃 GPU
- 如果同机同卡同时跑，容易 OOM
- 分阶段的好处：
  - 阶段边界清晰，失败可重跑某一段
  - 资源峰值可控（特别是 Boltz）
  - 更适合生产调度（夜间跑 Boltz，白天跑 Evo）

对应代码位置：
- `main.py` 的 `ArgumentParser(--stage)`
- 三个独立函数：`run_evo_stage`, `run_boltz_stage`, `run_analyze_stage`

### 5.2 `logic_type` 如何控制逻辑走向

`logic_type` 是业务规则路由器：

- `POSITIVE`
  - 生物学含义：某些突变可能使靶点“更能耐药”
  - 代码处理：用 `lethal_threshold` 决定是否淘汰

- `NEGATIVE`
  - 生物学含义：功能崩溃本身就可导致耐药
  - 代码处理：`delta < -8.0` 直接 `LOSS_OF_FUNCTION`，短路跳过 Boltz

- `STRUCTURAL`
  - 生物学含义：通道/结构完整性导向
  - 代码处理：双阈值区间区分致死 vs 功能坍塌

实现位置：
- `core/analyzer.py` 的 `apply_differential_filtering`

### 5.3 非编码区与可配置窗口

除了 CDS 区域，系统支持通过 `config.yaml` 配置非编码区窗口（`upstream`, `downstream`, `around_start`, `intergenic`）。
- **生成阶段**：`generator.py` 会根据链方向（`+` 或 `-`）正确截取序列，并生成突变。
- **Evo 阶段**：非编码区突变使用 `source_sequence` 作为野生型基准进行评分。
- **Boltz 阶段**：非编码区突变（`region_type='REGULATORY'`）会直接跳过结构模拟，因为它们不产生氨基酸变化。

### 5.4 多基因协同引擎 (Synergy Engine)

在 `analyze` 阶段，系统不仅评估单基因突变，还会通过 `PathwayAnalyzer.apply_multi_gene_synergy` 评估多基因协同效应。
- **规则驱动**：在 `config.yaml` 的 `synergy_rules` 中定义规则（如 `CO_LOF`, `CO_COMPLETED`, `REGULATOR_EFFECTOR_DEREPRESSION`）。
- **批次评估**：引擎会收集同一批次（或样本）中所有基因的突变状态，匹配规则。
- **得分修正**：匹配成功后，会为相关突变附加 `synergy_score` 和 `synergy_tags`，最终影响耐药性判定。

当前规则字段兼容两套命名（便于平滑升级）：
- 规则类型：`logic` 或 `rule_type`
- 分数修饰：`score` 或 `synergy_score_modifier`
- 解释文本：`interpretation` 或 `description`
- 效应子：`effectors`（列表）或 `effector`（单值）

引擎还支持：
- `tags`：命中后附加到 `synergy_tags`
- `apply_status` + `eligible_current_statuses`：命中后按条件覆盖状态

---

## 6. 关键调用链（方便你断点调试）

### 初始化链
- `scripts/init_project.py: main()`
  - `extract_target_gene_sequences`
  - `upsert_gene`
  - `MutationEngine.generate_saturation_mutagenesis`
  - `MutationEngine.bulk_insert_mutations`

### Evo 链
- `main.py: run_evo_stage`
  - `Evo2Client.get_likelihood`
  - `apply_differential_filtering`

### Boltz 链
- `main.py: run_boltz_stage`
  - `Boltz2Client.predict`
  - `_build_payload` → 三类 payload builder

### Analyze 链
- `main.py: run_analyze_stage`
  - `PathwayAnalyzer.apply_epistasis_corrections`
  - `PathwayAnalyzer.apply_multi_gene_synergy`

---

## 7. 新手读源码顺序（建议）

建议按下面顺序读，理解成本最低：

1. `config/config.yaml`（先懂业务数据）
2. `main.py`（懂总流程）
3. `core/analyzer.py`（懂判定规则）
4. `core/clients.py`（懂模型调用）
5. `scripts/init_project.py`（懂初始化）
6. `core/generator.py`（懂突变生成细节）
7. `core/database.py`（懂约束和迁移）

---

## 8. 你现在可以做的三件事（掌控项目）

1) 先运行一遍最小流程，观察每个阶段状态数量变化：
```sql
SELECT status, COUNT(*) FROM mutations GROUP BY status;
```

2) 随机抽一条 atpE 突变，从 `PENDING` 一路追到 CSV，验证“一个突变的一生”。

3) 改一个阈值（如 `collapse_threshold`），重跑 `--stage evo`，看分流变化，建立直觉。

---

## 9. 容器环境管理和路径自动探测（env_manager.py）

在全容器化部署中，`core/env_manager.py` 负责管理跨平台路径差异。

### 9.1 为什么需要环境管理

在容器和宿主系统中，数据路径完全不同：

**容器环境**：
- 本地数据：`/app/data`
- 共享存储：`/app/share_data`
- 日志：`/app/logs`

**宿主环境**（本地 Python 运行）：
- 本地数据：`/root/amr_hunter/data`
- 共享存储：`/root/gpufree-share/amr_hunter`
- 日志：`/root/amr_hunter/logs`

`env_manager.py` 的任务就是：
1. 自动检测运行环境（容器 vs 宿主）
2. 根据环境选择正确的默认路径
3. 智能探测共享存储位置
4. 提供统一的路径 API（不管你在容器还是宿主，调用方式一样）

### 9.2 自动探测共享路径的优先级

```python
# 核心逻辑（简化版）
def _resolve_share_path_container(self) -> Path:
    # 第 1 层：用户显式指定
    if env_share := os.getenv(self.SHARE_PATH_VAR):
        return Path(env_share)
    
    # 第 2 层：逐一检查候选路径是否存在
    for candidate in ["/app/share_data", "/mnt/shared", "/data/shared"]:
        if Path(candidate).exists():
            return Path(candidate)
    
    # 第 3 层：都没找到，使用默认值
    return self.CONTAINER_SHARE_ROOT
```

**容器环境候选顺序**：
1. `$AMR_SHARE_PATH`（如果设置）
2. `/app/share_data`（Docker 标准挂载）
3. `/mnt/shared`（网络存储常见路径）
4. `/data/shared`（数据目录备选）
5. `/app/share_data`（默认值）

**宿主环境候选顺序**：
1. `$AMR_SHARE_PATH`（如果设置）
2. `/root/gpufree-share/amr_hunter`（项目标准）
3. `/mnt/amr_share`（网络存储）
4. `/data/amr_share`（数据目录备选）
5. `./shared_data`（相对路径，同工作目录）
6. `/root/gpufree-share/amr_hunter`（默认值）

### 9.3 使用 EnvManager API

在任何模块中：

```python
from core.env_manager import get_env_manager

# 获取全局实例（单例）
manager = get_env_manager()

# 获取关键路径
db_path = manager.get_db_path()           # → /app/data/db/amr_mutations.db
results_dir = manager.get_results_dir()   # → /app/data/results
share_dir = manager.get_share_path()      # → 自动探测/默认的共享路径

# 跨容器同步
manager.sync_to_share(source_path, subdir="results")

# 自愈降级
fallback_path = manager.fallback_to_local(share_path)  # 共享存储不可用时本地化

# 调试：打印当前环境信息
manager.print_environment_info()
```

### 9.4 容器环境检测机制

```python
def is_running_in_container() -> bool:
    """检测是否在 Docker 容器内"""
    return os.path.exists("/.dockerenv")  # 容器内必有此文件
```

如果环境没有 `/.dockerenv`，也可显式指定：
```bash
export AMR_PLATFORM_TYPE=container
python main.py --stage all
```

### 9.5 修改自动探测候选路径

如果项目部署的标准路径变了，编辑 `core/env_manager.py`：

```python
CANDIDATE_SHARE_PATHS_CONTAINER = [
    Path("/app/share_data"),
    Path("/mnt/my_custom_nfs"),    # ← 添加你的路径
    Path("/data/shared"),
]
```

然后重建容器：
```bash
docker compose build amr-app
```

### 9.6 跨容器数据同步流程（auto_run.py 中的实际使用）

1. **启动时**（从共享存储恢复）：
   ```python
   DataSync.sync_share_to_local()  # /app/share_data → /app/data
   ```

2. **执行完毕**（同步回共享存储）：
   ```python
   DataSync.sync_local_to_share()  # /app/data → /app/share_data
   ```

这样即使容器停止，下次重启也能从共享存储恢复进度。

---

## 10. 维护注意事项（非常重要）

- 改 `logic_type` 或 `scenario` 时，先检查：
  - `config.yaml` 是否有值
  - `database.py` 的 CHECK 约束是否允许
  - `clients.py` payload builder 是否支持

- 改状态机时，必须同时检查：
  - `database.py` 状态约束
  - `main.py` 各阶段 SQL WHERE 条件
  - `analyzer.py` 的规则判断

- 改 `synergy_rules` 时，确保字段名与引擎兼容：
  - 推荐新写法：`rule_type`, `synergy_score_modifier`, `description`, `effectors`
  - 兼容旧写法：`logic`, `score`, `interpretation`, `effector`

- 不要直接改生产库结构，优先通过 `DatabaseManager._migrate_schema_if_needed` 追加迁移逻辑。

- **容器部署时的路径检查**：
  - 如改了 Docker 卷挂载点，务必同步更新 `CANDIDATE_SHARE_PATHS_*` 列表
  - 用 `docker compose logs amr-app` 观察路径探测日志，确保找到预期路径
  - 跨主机共享存储时，确保 `AMR_SHARE_PATH` 指向网络挂载（NFS/SMB），否则数据丢失
  - 本地模式（宿主 Python）测试时，验证 `EnvManager.env_name == "host"`

---

**相关文档**：
- [README.md](README.md) - 用户快速开始和故障排查
- 本文档 - 源码原理和架构理解
