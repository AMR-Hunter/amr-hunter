# 判定规则矩阵

## 范围

本文档记录当前每个基因、每类区域的表型判定覆盖情况。
这里把状态分成三类：

- 已实现：已经直接影响最终表型输出。
- 仅辅助：原始证据会保留，但不会直接控制最终表型。
- 暂不支持：当前流程可以存储或排序这类候选，但证据还没有校准到可稳定使用的规则层。

重要说明：本矩阵描述的是“最终表型判定支持度”，不是“候选是否能被流程处理”。
很多非编码候选已经可以被注释、入库、在 Evo 阶段保留，甚至在少数场景下进入 Boltz 排队；但如果缺少场景匹配的机制模型或校准阈值，仍然会在这里被标记为“暂不支持”或“仅辅助”。
另外，受支持的 regulatory PROTEIN_DNA 路径在默认配置下也没有全量打开运行，而是受 selection.include_regulatory 控制。

## 全局判定策略

| 领域 | 当前策略 | 状态 | 下一步 |
| --- | --- | --- | --- |
| 正向靶点 CDS 判定 | 需要 Evo 结果加上 Boltz 结构可行性；目前只有 atpE 有亲和力感知的逃逸逻辑 | 部分实现 | 用基因特异规则替代通用兜底逻辑 |
| 负调控因子 CDS 判定 | LOF 或结构塌陷可推断耐药；Rv0678 另外新增了面向受支持 PROTEIN_DNA 行的第一版调控区 DNA 结合逃逸规则 | 已实现 | 用真实分布校准 PROTEIN_DNA 阈值 |
| 外排泵结构型 CDS 判定 | 结构塌陷意味着敏感性增强或再敏化背景；mmpL5 和 mmpS5 现在也有第一版 COMPLEX 界面破坏规则，但只有在存在核心证据时才触发 | 已实现 | 校准 COMPLEX 阈值，并决定是否扩展到 mtrAB |
| 调控区判定 | 只有 PROTEIN_DNA 被视为受支持的调控区场景；其他调控区组合都会被显式标记为暂不支持 | 部分实现 | 加入基因特异的启动子或操纵子规则，停止混用非 DNA 场景 |
| 基因间区判定 | 目前通过锚定基因上下文转发处理；没有独立的基因间区判定模型 | 暂不支持 | 增加显式的启动子、操纵子、基因间区规则家族 |
| 结合类 Boltz 输出 | `binding_affinity`、`affinity_pred_value`、`pair_energy`、`complex_energy` 目前都只作为辅助原始输出持久化 | 部分实现 | 只有在按场景完成校准后才提升为正式判定依据 |

## 基因 × 区域矩阵

| 基因 / 分组 | 区域 | 配置场景 | 当前判定状态 | 当前规则 | 建议 |
| --- | --- | --- | --- | --- | --- |
| atpE | CDS | PROTEIN_LIGAND | 已实现 | 只有当 Evo 损伤强、Boltz 结构仍可行，且相对 WT 的亲和力偏移超过阈值时才判耐药；缺失亲和力时返回 `ATPE_AFFINITY_REQUIRED` | 作为主要正向靶点规则保留，并回填历史亲和力 |
| atpE | 调控区 | PROTEIN_LIGAND | 暂不支持 | 由于调控区证据不是 PROTEIN_DNA 场景，因此显式标记为暂不支持 | 不要从启动子变异直接推断 target escape；后续单独增加转录调控规则 |
| atpE | 基因间区 | 锚定基因派生 | 暂不支持 | 没有专门的基因间区判定逻辑 | 如果有需要，再增加启动子或基因间区表达模型 |
| atpB | CDS | PROTEIN_LIGAND | 仅辅助 | 已明确降级为 `ATPB_AUXILIARY_ONLY`；目前没有直接耐药规则 | 需要先决定 atpB 是直接耐药决定基因，还是仅作为 ATP 合酶复合体背景 |
| atpB | 调控区 | PROTEIN_LIGAND | 暂不支持 | 作为非 DNA 调控场景被显式标记为暂不支持 | 与 atpE 调控区相同：要么建表达效应模型，要么维持暂不支持 |
| atpB | 基因间区 | 锚定基因派生 | 暂不支持 | 没有专门的基因间区规则 | 只有在存在真实启动子或操纵子机制时再加入 |
| Rv0678 | CDS | PROTEIN_DNA | 已实现 | NEGATIVE 逻辑加上 Boltz 核心塌陷可产生 LOF 驱动的耐药判定 | 保持为当前调控因子 CDS 规则 |
| Rv0678 | 调控区 | PROTEIN_DNA | 已实现（第一版） | 对受支持的调控区行，当 Boltz 复合体置信度降低但整体折叠置信度仍保留时，触发 `RV0678_DNA_BINDING_ESCAPE` | 等更多 PROTEIN_DNA 证据后重新校准阈值 |
| Rv0678-mmpS5 / Rv0678-mmpL5 基因间区 | 基因间区 | 锚定基因派生 | 暂不支持 | 还没有单独的基因间区规则家族 | 增加显式的操纵子或基因间区调控模型 |
| mtrA | CDS | FOLDING_ONLY | 已实现 | STRUCTURAL 语义：LoF/结构塌陷 → 敏化（MtrAB 必需双组分系统，CRISPRi 敲低 → BDQ 超敏）；由 NEGATIVE 迁移而来 | 保持；如需复合体界面规则再评估 |
| mtrA | 调控区 | FOLDING_ONLY | 暂不支持 | 调控区变异在 PROTEIN_DNA 之外目前都不支持 | 将调控控制和蛋白复合体控制分开建模 |
| mtrA-mtrB 基因间区 | 基因间区 | 锚定基因派生 | 暂不支持 | 没有基因间区判定模型 | 只有在生物学上有充分理由时再增加启动子或操纵子模型 |
| mtrB | CDS | FOLDING_ONLY | 已实现 | 与 mtrA 相同：STRUCTURAL 语义，LoF → 敏化 | 与 mtrA 相同 |
| mtrB | 调控区 | FOLDING_ONLY | 暂不支持 | PROTEIN_DNA 之外的调控区变异目前不支持 | 与 mtrA 调控区相同 |
| mmpS5 | CDS | COMPLEX | 已实现（证据门控） | STRUCTURAL 塌陷会参与外排泵塌陷推断；当存在 COMPLEX 核心证据时，低界面分数且整体高置信度的结果会触发 `EFFLUX_COMPLEX_INTERFACE_DISRUPTED` | 在拿到更多 mmpS5 COMPLEX 证据之前保持保守 |
| mmpS5 | 调控区 | COMPLEX | 暂不支持 | PROTEIN_DNA 之外的调控区变异目前不支持 | 只有在有证据时才增加启动子特异规则 |
| mmpL5 | CDS | COMPLEX | 已实现（证据门控） | 与 mmpS5 相同，并在存在核心证据时使用第一版 COMPLEX 界面破坏规则 | 在更大的 mmpL5 证据集上重新校准阈值 |
| mmpL5 | 调控区 | COMPLEX | 暂不支持 | PROTEIN_DNA 之外的调控区变异目前不支持 | 与 mmpS5 调控区相同 |
| Rv1979c | CDS | FOLDING_ONLY | 已实现（探索性） | 探索性基因：证据照常入库，判定输出 UNKNOWN + EXPLORATORY_ONLY（CRyPTIC 校准 LoF 中性，方向未定） | 保留探索性定位，方向明确后再评估正式规则 |
| Rv1979c | 调控区 | FOLDING_ONLY | 暂不支持 | PROTEIN_DNA 之外的调控区变异目前不支持 | 只有在建立调控机制模型后才加入 |
| pepQ | CDS | FOLDING_ONLY | 已实现 | 通用 NEGATIVE 逻辑或结构塌陷语义 | 保留为第一版仅折叠规则 |
| pepQ | 调控区 | FOLDING_ONLY | 暂不支持 | PROTEIN_DNA 之外的调控区变异目前不支持 | 只有在建立调控机制模型后才加入 |
| glpK | CDS | FOLDING_ONLY | 已实现（探索性） | 与 Rv1979c 相同：UNKNOWN + EXPLORATORY_ONLY | 与 Rv1979c 相同 |
| glpK | 调控区 | FOLDING_ONLY | 暂不支持 | PROTEIN_DNA 之外的调控区变异目前不支持 | 只有在建立调控机制模型后才加入 |

## 已解决事项

| 主题 | 当前状态 |
| --- | --- |
| atpE 缺失亲和力时的兜底逻辑 | 已在代码中解决：不再直接回退到耐药，而是要求必须有亲和力上下文 |
| 不支持的调控区场景 | 已在代码中解决：会显式标记为暂不支持，而不是悄悄落成普通 reviewed |
| atpB 的定位 | 已在代码中以第一版形式解决：明确只作为辅助证据，不作为已校准的直接耐药规则 |
| Rv0678 的 PROTEIN_DNA 规则家族 | 已以第一版形式落地：受支持的调控区行现在可以提升为 `RV0678_DNA_BINDING_ESCAPE` |
| mtrA/mtrB 判定方向修正 | 已解决：由 NEGATIVE（LoF → 耐药）迁移为 STRUCTURAL（LoF → 敏化）；`mtrAB_deregulation` 协同规则及判定分支已移除 |
| Rv1979c/glpK 探索性基因短路 | 已解决：判定引擎输出 `UNKNOWN + EXPLORATORY_ONLY`，不产出正式 R/S |
| 外排泵 COMPLEX 规则家族 | 已以第一版形式落地：针对 mmpS5 和 mmpL5，受支持的 CDS 行现在可以提升为 `EFFLUX_COMPLEX_INTERFACE_DISRUPTED` |

## 尚未解决事项

| 主题 | 仍然开放的原因 |
| --- | --- |
| PROTEIN_DNA 结合阈值 | Rv0678 的第一版阈值已经存在，但它们目前还是保守占位值，仍需用真实 PROTEIN_DNA 分布做校准 |
| COMPLEX 界面阈值 | mmpS5 和 mmpL5 的第一版外排泵阈值已经存在，但它们目前还是保守占位值，仍需用真实 COMPLEX 分布做校准 |
| 基因间区判定模型 | 基因间区目前仍依赖锚定基因上下文，还没有自己的机制感知规则 |
| atpB 的直接耐药模型 | 明确延后，直到选定基因特异的机制规则为止 |
| mtrAB 的 COMPLEX 界面规则 | 还没有实现；当前 COMPLEX 扩展优先先做在外排泵结构基因上 |

## 建议的下一轮迭代顺序

1. 回填历史 atpE protein-ligand 亲和力，并重新运行分析。
2. 重新运行或回填受支持的 Rv0678 PROTEIN_DNA 样本和受支持的 efflux COMPLEX 样本，让新的第一版规则真正作用到历史数据上。
3. 明确 atpB 是继续保持辅助定位，还是升级为直接承载规则的目标基因。
4. 在检查原始指标分布后，再决定是否把 COMPLEX 规则家族从外排泵结构基因扩展到 mtrA 和 mtrB。