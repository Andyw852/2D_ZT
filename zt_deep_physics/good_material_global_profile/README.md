# 高性能热电材料相对全局的必要特征分析

本目录回答一个比相关性更具体的问题：在 Starrydata2 的全体可用实验样品中，峰值 ZT 较高的材料通常落在什么物性范围，这些范围相对全局是否富集，以及哪些性质目前没有足够证据设定硬阈值。

## 分析口径

- 每个实验 `sample_id` 只保留一个峰值 ZT，其他性质只在同一样品的峰值 ZT 温度处插值配对；
- “高性能”主定义为全局峰值 ZT 的前 10%，同时检查前 25%、前 5% 和 `ZT >= 1`；
- `P10-P90` 是覆盖 80% 高性能样品的**典型窗口**，`P05-P95` 是覆盖 90% 的**宽筛选窗口**；
- 用分位数富集、Mann-Whitney AUC、KS 距离和完整样本覆盖率共同判断，不把单变量相关性写成因果关系；
- 不要求所有性质同时非空。每个性质都使用自己的最大成对完整样本集，并明确报告样本数。

这里的“必要”是经验意义的高覆盖条件，不是普适物理定律。单一性质通常也不是充分条件；真实设计仍需同时满足较高功率因子和较低总热导，并检查稳定性、可合成性与工作温区。

## 输入

默认读取：

`../empirical/outputs/experimental_ZT_with_structure_metadata.csv`

该表由 `../empirical/build_empirical_atlas.py` 从 Starrydata2 曲线和样品元数据构建。

## 运行

```bash
cd /home/wangchao/work_wc/2D_ZT
~/miniconda3/envs/te_manifold/bin/python \
  zt_deep_physics/good_material_global_profile/analyze_good_materials.py
```

## 产物

- `outputs/report.md`：中文结论、推荐窗口和证据边界；
- `outputs/global_vs_good_numeric.csv`：全局与高性能组的分位数、覆盖率和效应量；
- `outputs/screening_ranges.csv`：典型窗口、宽窗口和跨阈值稳健核心；
- `outputs/quantile_enrichment_bins.csv`：每个全局十分位区间的高性能富集倍数；
- `outputs/threshold_sensitivity.csv`：高性能定义变化时的范围稳定性；
- `outputs/joint_rule_performance.csv`：单性质与联合软规则的保留率、命中率和富集；
- `outputs/categorical_enrichment.csv`：样品形态和材料家族的富集结果；
- `outputs/mechanism_model_priors.csv`：项目深层物理模型的情景先验（明确不等同于实验全局范围）；
- `outputs/analysis_manifest.csv`：输入、阈值和样本数审计；
- `figures/01_good_ranges_on_global_percentiles.png`：高性能窗口在全局百分位中的位置；
- `figures/02_decile_enrichment_heatmap.png`：不同性质区间的高性能富集；
- `figures/03_effect_vs_coverage.png`：效应强度与数据覆盖率；
- `figures/04_threshold_robustness.png`：更换高性能定义后的区间稳定性。
