# matexplore —— 通用多智能体"假设→生成→验证→反馈"闭环发现框架

一个**领域无关**的材料/性质发现编排框架，第一个落地实例是**二维热电材料（2D thermoelectrics）**。

## 设计

```
        ┌──────────────────── 闭环(可多轮) ────────────────────┐
        ▼                                                      │
  [Analyst 分析]  读数据库 -> 计算"易算特征<->目标"相关 -> 提假设(折线图)
        │                                                      │
        ▼                                                      │
  [Generator 生成] 原型替换/组成枚举 -> 生成候选结构(POSCAR) + 替代模型打分
        │                                                      │
        ▼                                                      │
  [Validator 验证] T0 替代筛选(te-screen) -> T1 MACE 弛豫(3090) -> T2 DFT 输运(jzzn)
        │                                                      │
        ▼                                                      │
  [Collector 收集] 汇总结果 -> 对照假设 -> 反馈给 Analyst ──────┘
```

四个智能体各自独立、通过共享 `round_store` + `artifacts` 传递中间产物；
每个智能体的"推理"可插拔 LLM（`config.yaml` 的 `llm.provider`），缺省走确定性数据推理。

## 目录

```
matexplore/
├── config.yaml                  # 全局配置(领域相关的知识库路径/生成规则/超算)
├── requirements.txt
├── matexplore/                  # 框架包(领域无关骨架)
│   ├── orchestrator.py          #   闭环编排器
│   ├── config.py                #   配置加载(内置 YAML 子集解析，免 PyYAML 依赖)
│   ├── knowledge/               #   数据库读取 + 轮次持久化
│   ├── hypothesis/              #   特征计算 + 相关面板 + 折线图
│   ├── generation/              #   结构生成 + 替代模型打分
│   ├── validation/              #   taskflow(tf) 客户端
│   └── agents/                  #   analyst/generator/validator/collector
├── skills/te-screen/            #   taskflow 技能(同时安装到 ~/software/taskflow/skill/)
├── scripts/                     #   分析/训练/入口脚本
├── reports/                     #   相关面板/折线图/候选/模型指标
└── runs/                        #   每次闭环运行的结果(假设/候选/验证/反馈 + POSCAR)
```

## 快速开始

```bash
# 环境(需要 numpy/pandas/scipy/scikit-learn/matplotlib；本项目用 te_manifold)
conda activate te_manifold

# 1) 重建特征面板 + 相关折线图(纯本地分析)
python scripts/build_feature_target_analysis.py

# 2) (重)训练易算特征 -> ZT_e/PF 替代模型(输出到 skills/te-screen/)
python scripts/train_surrogate.py

# 3) 跑多轮多智能体闭环(默认 dry_run，不提交超算)
python scripts/run_pipeline.py --rounds 3

# 4) 对闭环胜者真实验证：提交 3090 MACE 弛豫(复用 opt-mace-gpu 引擎)
python scripts/run_3090_validation.py --n 5
```

## 验证三层（按开销递增）

| 层 | 技能 | 机器 | 算什么 | 用途 |
|---|---|---|---|---|
| T0 | `te-screen`(新开发) | 登录节点 | 易算特征 -> ZT_e/PF 预测 | 秒级粗筛，闭环内快速反馈 |
| T1 | `opt-mace-gpu` | 3090 GPU | MACE 弛豫 + 形成能 | 稳定性第一道闸 |
| T2 | `band-dft-cpu`+`ke-dft-cpu` | jzzn CPU | Eg/m*/S/σ/PF/κe (DFT) | 真值确认 |

> 晶格热导率 `kl-*`（κ_L）本轮按用户要求**跳过**（大体系暂不验证）；
> 它是"ZT_e 上限 → 真实 ZT"升级的直接路径，见复盘报告 L1。

## 复用到其它领域

框架的骨架（orchestrator / agents / round_store / 相关面板 / 折线图）与热电无关。
换领域只需：
1. 提供一份"材料/样本 → 易算特征 + 目标标签"的 `panel_csv`；
2. 在 `hypothesis/cheap_features.py` 实现该领域的特征提取（或直接给特征表）；
3. 在 `generation/` 换生成器（如分子生成模型 / 组合库枚举）；
4. 在 `validation/` 接上该领域的计算引擎（taskflow 技能或其它）。

## 结构生成模型说明

本环境未安装 MatterGen / CDVAE / DiffCSP / FlowMM 等生成式晶体结构模型，
故 `generation/structure_models.py` 提供可靠的离线实现：**原型替换（元素族同价替换）+ 组成枚举**。
生成式模型可通过 `generate_from_model()` 接口接入（见文件内注释）。

## 多轮闭环（已实现）

`runs/<run_id>/ledger.json` 是跨轮累积记忆，让闭环真正闭合：
- **去重**：生成器跳过已生成过的公式（`ledger.seen`）；
- **围绕胜者迭代**：collector 把本轮 T0 胜者（含结构）写回 `ledger.winners`，下一轮生成器以「基线种子 + 上一轮胜者种子」继续替换，候选逐步精炼；
- **验证状态累积**：`ledger.validated` 记录每个候选 T0/T1 状态。

实测 3 轮：130 个唯一公式，每轮 top 候选随胜者前移而演化（F4Sn → ClF3Pb → Cl2F2PbSn 等）。

## 3090 真实验证（已跑通）

`scripts/run_3090_validation.py` 忠实复用 taskflow 技能 `opt-mace-gpu` 的引擎 `mace_relax.py` +
模型 `MACE-matpes-pbe-omat-ft`，ssh 到 3090 以 mace-gpu 环境 + CUDA 弛豫，拉回 `relax_summary.json`。
幂等：已 `T1_mace_relaxed` 且 converged 的候选自动跳过。

实测 top-5 全部收敛（力 <1e-3 eV/Å、应力 <0.002 GPa），E/atom 排序：
ClF3Pb2(-4.026) < Cl2F2PbSn(-3.883) < Cl2F2GePb(-3.882) < BrF3Pb2(-3.874) < BrClF2PbSn(-3.787) eV。
结果见 `runs/<run_id>/validation_3090_all.csv` 并回写 ledger。

> 注意：MACE 弛豫 converged 只说明"在 MACE 势里存在局域极小"（弱稳定性信号），
> 不等于热力学稳定（需凸包/形成能）也不等于动力学稳定（需声子谱无虚频）。

## 诚实边界

- 替代模型（易算特征→ZT_e）CV Spearman ≈ +0.35(n)/+0.49(p)：只适合**粗筛/排名**，不是精确预测；
- 未含 κ_L，故所有"ZT"都是电子上限 ZT_e；
- 候选的稳定性（E_hull 全 0）未经 DFT 确认，必须经 T1 稳定性闸过滤。
