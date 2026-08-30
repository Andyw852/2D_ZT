# 二维热电优值的深层物理量拆解与结构设计

这个目录是对根目录现有 `charts/` 分析的一次独立重构。旧文件不修改。项目现在分为两层：

1. **主结果：`empirical/` 总体数据图。** 尽量使用 Starrydata2、JARVIS 2D 和已有跨库匹配中的全部可用样品，每个面板按其所需字段独立取最大完整样本集；
2. **辅助结果：原有 `figures/` 参数模型图。** 使用自洽的二维单抛物带模型、可变 Lorenz 数和显式声子模态近似解释总体统计，但不代表材料总体分布。

新增的 `good_material_global_profile/` 在主结果之上比较高 ZT 样品与全局分布，给出软筛选范围、区间富集、阈值敏感性和证据覆盖率，并把实验范围与参数模型先验严格分开。

两层共同回答以下问题：

1. `ZT = S²σT/(κe+κL)` 如何继续拆成更深的可计算物理量；
2. 载流子浓度、两类有效质量、群速度和声子寿命的可行数值窗口；
3. 褶皱、孔隙/骨架、软单元、内部成键和声—光支关系如何进入这些参数；
4. 哪些结构描述可以定量，哪些应删掉或改写；
5. 二维响应面与局部二次模型是否能给出稳定的设计区间。

> 当前 README 先记录模型和复现入口。运行结果和图表说明由 `src/run_analysis.py` 自动更新到 `outputs/summary.md`。

## 目录

```text
zt_deep_physics/
├── config/baseline.json              # 基准情景与扫描范围
├── data/feature_definitions.csv       # 结构特征的保留、替换、删除依据
├── figures/                           # 生成的 5 张图
├── empirical/                         # 主结果：全部可用材料的总体统计图
├── good_material_global_profile/      # 高 ZT 组相对全局的范围、富集和筛选规则
├── outputs/                           # 数值表、最优窗口、二次拟合系数、结论
├── src/zt_model.py                    # 电子与晶格物理模型
├── src/run_analysis.py                # 扫描、拟合、绘图、报告
└── tests/test_zt_model.py             # 量纲、极限和单调性测试
```

## 核心分解

二维电子侧采用未归一化 Fermi–Dirac 积分

\[
F_j(\eta)=\int_0^\infty\frac{x^j}{1+\exp(x-\eta)}\,dx.
\]

对二维抛物带和声学形变势散射（`r=0`），输运谱指数 `a=r+d/2=1`：

\[
n_{2D}=\frac{N_v m_d^* k_BT}{\pi\hbar^2}F_0(\eta),
\]

\[
|S|=\frac{k_B}{e}\left[\frac{2F_1}{F_0}-\eta\right],\qquad
L=\left(\frac{k_B}{e}\right)^2\left[\frac{3F_2}{F_0}-\left(\frac{2F_1}{F_0}\right)^2\right],
\]

\[
\mu_{2D}=\frac{e\hbar^3 C_{2D}}{k_BT\,m_c^*m_d^*E_1^2},\qquad
\sigma=\frac{n_{2D}}{t_{eff}}e\mu,\qquad \kappa_e=L\sigma T.
\]

这里 `m_d*` 是每个能谷的态密度质量，`m_c*` 是指定输运方向的导电质量，`Nv` 是谷简并度。把二者混成一个“有效质量”会隐藏最重要的带结构设计自由度。

晶格侧用模态 BTE 的压缩形式：

\[
\kappa_{L,\alpha\alpha}=\frac{1}{V}\sum_\lambda C_\lambda v_{\lambda,\alpha}^2\tau_\lambda
\approx \frac{C_V v_{g,\alpha}^2\tau_{eff}}{2}f_{pore}.
\]

其中二维面内各向同性平均给出分母 2；寿命用 Matthiessen 规则合并本征、孔边界、褶皱和声—光支重叠散射。该压缩模型用于趋势与范围推演，不替代 DFPT + ShengBTE/Phono3py。

## 复现

项目使用已有 `te_manifold` Conda 环境：

```bash
cd /home/wangchao/work_wc/2D_ZT
~/miniconda3/envs/te_manifold/bin/python zt_deep_physics/src/run_analysis.py
~/miniconda3/envs/te_manifold/bin/python zt_deep_physics/empirical/build_empirical_atlas.py
~/miniconda3/envs/te_manifold/bin/python zt_deep_physics/good_material_global_profile/analyze_good_materials.py
~/miniconda3/envs/te_manifold/bin/python -m pytest -q \
  zt_deep_physics/tests zt_deep_physics/good_material_global_profile/test_analysis.py
```

## 证据边界

- `精确关系`：ZT 恒等式、Fermi 积分 SPB、`κe=LσT`、模态 BTE。
- `低阶近似`：二维形变势迁移率、有效厚度换算、压缩后的单一 `Cv-vg-τ` 晶格模型。
- `情景参数化`：孔隙、褶皱、软单元和声—光支重叠的散射强度。输出只能解释方向、耦合和条件窗口，不能冒充具体材料预测。
- 二维材料的 `W m⁻¹ K⁻¹` 和 `S m⁻¹` 都依赖厚度定义；因此同时输出片密度 `cm⁻²` 和体等效浓度 `cm⁻³`。

## 主要文献锚点

- Zhu et al., [Restructured single parabolic band model for quick analysis in thermoelectricity](https://doi.org/10.1038/s41524-021-00587-5), *npj Computational Materials* 7, 116 (2021).
- Qiao et al., [High-mobility transport anisotropy and linear dichroism in few-layer black phosphorus](https://doi.org/10.1038/ncomms5475), *Nature Communications* 5, 4475 (2014).
- Li et al., [ShengBTE: A solver of the Boltzmann transport equation for phonons](https://doi.org/10.1016/j.cpc.2014.02.015), *Computer Physics Communications* 185, 1747 (2014).
- Lee et al., [Resonant bonding leads to low lattice thermal conductivity](https://doi.org/10.1038/ncomms4525), *Nature Communications* 5, 3525 (2014).
- Feng et al., [Thermal Conductivity of Graphene Wrinkles](https://doi.org/10.1021/acs.jpcc.6b07162), *J. Phys. Chem. C* 120, 23807 (2016).
