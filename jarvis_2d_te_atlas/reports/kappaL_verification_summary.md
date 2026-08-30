# kL (晶格热导率) 与结构 / 电子性质关联性验证

> **更正提示（2026-08-27）**：本报告初版的子集跨视图近邻比较使用了局部行号，且成对距离相关的显著性没有处理距离对之间的非独立性。两项问题已用材料 JID 近邻和节点置换（Mantel 型）检验修正；后续结论请以 [kappaL_dual_manifold_reanalysis.md](kappaL_dual_manifold_reanalysis.md) 为准。

> 复用既有「多视图相似性 + 跨视图近邻重叠」流程，回答：结构相似、电子性质相似，与晶格热导率 kL 相似之间是否存在直接联系。

## 1. 数据源核查结论

| 数据 | 位置 | 是否含 kL | 说明 |
|---|---|---|---|
| JARVIS dft_2d | 本地 data/raw/jarvis | 否 | 仅电子输运(Seebeck/PF/cond/电子 ke)，无晶格 kL |
| JARVIS dft_3d | 本次下载(figshare jdft_3d-8-18-2021.json, 55723 条) | 否 | 同 dft_2d，nkappa/pkappa 为电子 ke；全字段扫描无 kl/Slack/声子项 |
| MP Phonon v1.1 (Zenodo) | 未采用 | 否 | 原始 phonopy 数据，无直接 kL 表 |
| Materials Project Slack kL | 需 API key | 待用 | 最干净的 DFT kL，但当前环境无 key (401) |
| AFLOW kL | API | 不可用 | 服务 503 / TLS 异常 |
| starrydata2 (本地已下载) | data/raw/external/starrydata2 | 是 | 实验 Lattice thermal conductivity 曲线 6740 条 |

结论：本地没有「结构 + 电子性质 + DFT 晶格热导率」三合一现成库。采用 实验 kL(starrydata2, 本地) + DFT 结构/电子(JARVIS dft_3d, 本次下载) 的跨库组合，
按 pymatgen 归一化化学式精确映射，得到 137 个共有稳定材料(每化学式选最稳定多晶型)。

## 2. 视图构造 (与既有 2D 流程一致)

- Structure：几何-only SOAP(dummy 物种 X, n_max=6/l_max=6, mean 池化, r_cut=6) + 组成分数 Hellinger 距离, d = 0.5*d_geo + 0.5*d_comp (scripts 20/22 同款)。
- Electronic：OptB88vdW 带隙 + 电子/空穴有效质量 (RobustScaler)。
- kL：starrydata2 在 300K 的晶格热导率(同成分取中位数)，取 log10。

## 3. 结果一：跨视图近邻重叠 (k=10, 1000 次置换基线)

| 视图对 | N | 重叠 | 随机基线 | z | p |
|---|---|---|---|---|---|
| Structure 对 kL | 137 | 0.105 | 0.073 | 4.2 | <1e-3 |
| Structure(几何) 对 kL | 137 | 0.102 | 0.072 | 3.9 | <1e-3 |
| Structure(组成) 对 kL | 137 | 0.103 | 0.073 | 4.2 | <1e-3 |
| Eg 对 kL | 137 | 0.074 | 0.076 | -0.4 | 0.63 |
| Electronic 对 kL | 82 | 0.082 | 0.071 | 1.4 | 0.098 |
| Elastic(B,G) 对 kL | 82 | 0.082 | 0.072 | 1.1 | 0.15 |

距离相关 (Spearman, 成对距离)

| 视图对 | Spearman |
|---|---|
| Structure 对 kL | +0.220 |
| 几何-only 对 kL | +0.195 |
| 组成-only 对 kL | +0.082 |
| Eg 对 kL | +0.095 |
| Electronic 对 kL | +0.075 |
| Elastic(B,G) 对 kL | +0.199 |

## 4. 结果二：kL 与单个物理描述符的直接相关 (log10 kL 的 Spearman)

| 描述符 | N | Spearman | 解读 |
|---|---|---|---|
| Eg (带隙) | 137 | -0.011 | 无关 |
| m_elec / m_hole | 82 | -0.02 / -0.04 | 无关 |
| v_s 代理 = sqrt(B/rho) | 82 | +0.604 | 强正相关(声速, Slack 机制) |
| B_kv (体模量) | 82 | +0.529 | 强正相关 |
| G_gv (剪切模量) | 82 | +0.505 | 强正相关 |
| avg_mass (平均原子质量) | 137 | -0.403 | 越重 kL 越低 |
| density | 137 | -0.149 | 弱负相关 |

## 5. 结论

1. 电子性质 -> kL：无直接联系。带隙、有效质量与 kL 直接相关约 0 (Spearman -0.01~-0.04)，
   近邻重叠落在随机基线(Eg 对 kL 的 z=-0.4)。与此前「电子-电子输运共享公共空间」正交，物理自洽：
   晶格热导是声子输运，与电子能带结构无关。

2. 结构/弹性 -> kL：存在真实但集中于特定通道的联系。结构相似性对 kL 近邻有弱但显著的正向保持
   (重叠 z约4、成对距离 Spearman +0.22)，根本驱动是弹性/声速/原子质量这条 Slack 通道：
   sqrt(B/rho) 与 kL 的 Spearman 达 +0.60，体模量 +0.53、剪切模量 +0.51、平均原子质量 -0.40。
   即不是广义结构相似直接等于 kL 相似，而是结构中决定声子群速度的因素(键刚性、质量)在起作用。

3. 与此前结论的呼应：电子输运 <-> 电子能带(公共空间)；晶格输运 <-> 结构/弹性(本验证)。
   两条输运通道的解耦在此得到独立印证：kL 与电子视图无关、与结构/弹性视图相关，
   恰是电子输运结论的镜像。

## 6. 局限与后续

- 样本 137(有效质量/弹性子集 82)，且 kL 为实验值(含掺杂、合成、测量噪声，取同成分中位数)，
  同成分对映最稳定多晶型结构有一定歧义。
- 升级路径：若提供 Materials Project API key，可用其 Slack 模型 kL(~5k 材料) 做纯 DFT、大样本复验
  (结构 + Eg + 有效质量 + kL 全 DFT 同源)，预计 structure/弹性->kL 更稳、electronic->kL 仍约 0。

## 产物

- 特征/视图: features/kl_verify/kl_views.parquet, data/processed/kl_soap_geo.npy, kl_comp_frac.npy
- 审计表: data/audit/kl_view_overlap.csv, kl_view_distance_corr.csv, kl_descriptor_corr.csv
- 图: figures/kl_descriptor_scatter.png, figures/kl_struct_vs_kL_dist.png
- 脚本: scripts/kl_verify_01_recon.py ~ 04_descriptors.py
- 原始数据: data/raw/external/jarvis_kl/jdft_3d-8-18-2021.json(.zip)
