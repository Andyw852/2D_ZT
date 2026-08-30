# 全体材料数据图（主分析）

这一层按“每个关系使用它所需的最大完整样本集”构图，不要求所有样品同时具备 ZT、载流子浓度、Seebeck、电导和热导。

## 配对原则

- 真实 ZT 总体图：每个 Starrydata2 `sample_id` 只保留一个有效峰值 ZT；其他性质只从同一 `sample_id` 插值到峰值温度。
- 不要求 ZT 的电子/热输运图：分别在 300、600、900 K 对同一样品的两种性质作温度内插配对。
- 晶格热导温度带：每个温度先对同一样品的重复曲线取中位数，再统计样品分位数，避免点数多的曲线获得额外权重。
- 结构元数据：相对密度、晶粒尺寸和样品形态来自 Starrydata2 样品元数据，和真实 ZT 按 `sample_id` 合并。
- 二维形状代理：JARVIS 2D 本身没有真实 ZT 和孔隙/褶皱标签；使用全部可用固定条件 Seebeck/电导，明确标为 DFT 数据而非实验 ZT。
- 声速代理：实验 κL 与 JARVIS 3D 弹性数据按化学式跨库匹配，存在晶型和实验状态歧义，因此单独标记。

## 产物

- `figures/01_data_coverage.png`：数据覆盖和各关系配对数；
- `figures/02_experimental_ZT_global.png`：真实 ZT 与 n、S、σ、κL 的总体分布；
- `figures/03_electronic_all_available.png`：无需 ZT 的全部电子输运配对；
- `figures/04_thermal_all_available.png`：无需 ZT 的全部热输运配对；
- `figures/05_structure_metadata_experimental.png`：相对密度、孔隙、晶粒和样品形态；
- `figures/06_jarvis2d_shape_transport.png`：二维原子几何代理与 S/σ；
- `outputs/panel_manifest.csv`：每个面板的数据源、样本数、温度和配对规则；
- `outputs/empirical_summary.md`：中文结论与限制。
- `FIGURE_GUIDE.md`：六张图、20 个子图的中文逐图解释、读图规则和证据边界。

运行：

```bash
cd /home/wangchao/work_wc/2D_ZT
~/miniconda3/envs/te_manifold/bin/python zt_deep_physics/empirical/build_empirical_atlas.py
```
