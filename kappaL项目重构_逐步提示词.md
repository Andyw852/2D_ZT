# kappaL-multiview-verify 审计与重构：逐步提示词文档

> 用法：每一步是一个独立的 prompt，按顺序贴给 Claude Code（或你用的 agent），一步跑完、验收通过再进下一步。
> 不要一次把整份文档贴进去 —— 分步执行才能在每一步卡住错误。
> 文中 `[[ ]]` 包住的是你需要替换的内容。

---

## 目录

- 第 0 部分：立刻要做的两件事
- 第 1 部分：现状诊断（证据清单，供你判断要不要重构）
- 第 2 部分：重构后的科学问题
- 第 3 部分：Step 0 – Step 13 提示词
- 第 4 部分：数据源清单
- 第 5 部分：重构后 README 应该长什么样

---

# 第 0 部分：立刻要做的两件事

## 0.1 吊销泄露的 API key

`mp_kappaL/download_mp.py` 第 5 行明文提交了 Materials Project 的 API key，已经在 git 历史里。

现在就做：

1. 登录 Materials Project → Dashboard → 重新生成 API key（老 key 立即失效）
2. 新 key 只放环境变量，永远不写进代码
3. 代码里改成 `KEY = os.environ["MP_API_KEY"]`，取不到就直接报错退出

> 注意：改代码不会从 git 历史里删掉旧 key。历史清洗（`git filter-repo` / BFG）可以做，但既然 key 已经吊销，优先级不高。不要因为"想清洗历史"而拖延吊销。

## 0.2 给现状打个 tag 存档

```bash
git tag -a v0-original -m "重构前存档：原始多视图分析"
git push origin v0-original
```

这样重构过程中随时能对比，也保留了你已经做过的工作痕迹。

---

# 第 1 部分：现状诊断（证据清单）

下面每一条都是我在你 clone 下来的仓库和 `processed/` 数据上实测的，不是推测。你可以自己复现。

## 问题 1 —— 主结论是循环论证【致命】

**证据**（在 `views_meta.parquet` 上，先剔除物理不可能值后 N=12190）：

| 回归 | R² |
|---|---|
| `log10(clarke) ~ log(debye) + log(density)` | **0.979** |
| `log10(clarke) ~ log(v_long) + log(v_trans) + log(cahill)` | **0.994** |

Clarke 最小热导率的定义式是 κ_min = 0.87·k_B·M_a^(−2/3)·E^(1/2)·ρ^(1/6)，等价于 κ ∝ n^(2/3)·v_m；
Debye 温度的定义式是 Θ_D ∝ n^(1/3)·v_m。二者共享同一组变量，`Spearman(Debye, log κ) = +0.99` 是代数结果。

**后果**：README 第 1 节"弹性(力学) ↔ kL：强联系 +0.65，富集 25x"和第 4 节"Debye +0.99、纵声速 +0.96"都不是物理发现。
`reports/02` 里在 Debye 那一行小字标了"一致性检查"，说明你其实意识到了 —— 但 README 的头条表没有反映这个认识，读者会当成结论。

**改法**：这一条只能作为**代码正确性检查**（"我的流水线能否恢复已知的解析关系"），必须从结论表里移出去。

## 问题 2 —— 目标变量不是晶格热导率【致命】

`clarke` / `cahill` 是**非晶极限的最小热导率**，不是晶态 κ_L。

**证据**（清洗后 N=12190）：

| 分位 | p1 | p25 | p50 | p75 | p99 |
|---|---|---|---|---|---|
| clarke (W/mK) | 0.098 | 0.418 | 0.664 | 0.994 | 2.92 |

四分位区间只跨 **2.4 倍**。真实 κ_L 在这批材料里应跨 4 个数量级（金刚石 ~2000，笼状化合物 ~0.3）。

**后果**：你不是在研究"什么决定 κ_L"，而是在研究"什么决定 κ_min 这个公式的取值"。而 κ_min 公式的取值由声速和数密度决定 —— 回到问题 1。

**改法**：换真实 κ_L 数据源（见第 4 部分）。κ_min 只能作为参考下界。

## 问题 3 —— 度量没有标定，无法支撑任何"弱"或"无关"的判断【致命】

Elastic 与 kL_clarke 在数学上是 R²=0.98 的确定性关系，但你的跨视图 Spearman 只给出 **+0.65**。

**后果**：确定性关系的度量上限是 0.65 —— 那么 Structure 的 0.36 是"中等"还是"其实很强"？Eg 的 0.09 是"无关"还是"度量太钝"？**目前无法判断。**

原因之一是维度不对等：kL 视图是 1 维距离 |Δlog κ|，Elastic 是 3 维欧氏距离，Structure 是几百维。
高维距离矩阵存在测度集中（所有距离趋同），跨维度比较 Spearman 值不是同一把尺子。

**改法**：做**标定曲线** —— 人工构造已知依赖强度（R² = 0.1/0.3/0.5/0.7/0.9）的合成目标变量，跑同一套跨视图流程，画出"真实依赖强度 → 度量读数"的映射。有了这条曲线，0.36 才能被翻译成有意义的话。

## 问题 4 —— 电子视图的匹配严重污染【严重】

**证据**：

- `gap_opt` 中 **69.65%（7687/11037）恰好等于 0** → 该视图的距离矩阵有一个巨大的全零块，kNN 完全由 argpartition 的排序偶然性决定
- `graph_utils.kNN_affinity` 里专门写了 `tiebreak_seed` 处理这个问题，但 `crossview_analysis.py` 用的是自己定义的 `knn_sets`，**没有用**
- **14.8% 的 MP 材料共享同一个 JARVIS canon 化学式**：26 个 SiO₂ 多晶型、21 个 Al₂O₃、18 个 C、15 个 Si、15 个 SiC 各自被赋予**完全相同**的一组电子性质
- `download_mp.py` 已经下载了 MP 自己的 `band_gap`（按 material_id 严格对应，零匹配噪声）到 `summary_bandgap.json`，**后续从未使用**

**后果**："Eg ↔ kL 富集 3.4x"极可能是两个人造效应的叠加：(a) 重复点抬高重叠，(b) 金属/非金属二分与 κ 相关。跟"带隙相似"没关系。
而"电子性质与 kL 无关"这个结论是从**噪声极大的匹配**上得到的**零结果**，零结果在这种噪声水平下不可解释。

**改法**：电子性质必须按 material_id 对应，不能按化学式跨库匹配。用 MP 自己的 band_gap，或用 Ricci BoltzTraP 数据库（也是 MP-id 索引）。

## 问题 5 —— 脏数据没在源头清洗【严重】

**证据**（`views_meta.parquet`，N=12246）：

| 异常 | 条数 | 最大值 |
|---|---|---|
| `debye > 5000 K` | 46 | 3.48×10⁸ K |
| `shear_vrh > 1000 GPa` | 57 | 5.6×10¹³ GPa |
| `v_long > 2×10⁴ m/s` | 52 | 3.29×10⁹ m/s（**超光速**）|
| `bulk_vrh > 1000 GPa` | 46 | 1.37×10⁸ GPa |
| `snyder_acoustic > 10⁴` | 51 | 2.4×10¹⁹ |

`bulk_vrh < 1000` 的过滤只写在 `elastic_filt` 里（只作用于 Elastic 视图），Debye / 声速 / Structure↔kL / Eg↔kL 的计算全程带着这些值。

**改法**：清洗必须在 `build_views.py` 里做一次，写成显式的物理合理性区间，并输出被剔除样本的清单（可审计）。不要在下游脚本里各清各的。

## 问题 6 —— 三个可复现的代码 bug【严重】

**bug A：仓库开箱即坏。**
`crossview_analysis.py:14` 和 `descriptor_analysis.py:71`：
```python
sys.path.insert(0, str(root / "jarvis_2d_te_atlas" / "scripts"))
```
`jarvis_2d_te_atlas/` 不在仓库里，`graph_utils.py` 实际在 `scripts/`。clone 下来直接 ImportError。

**bug B：`descriptor_analysis.py` 行错位，第二张图的数值是错的。**
```python
meta = meta[meta["bulk_vrh"] < 1000].copy()   # 掉了 46 行
...
d_geo = soap_distance(soap_geo)               # soap_geo 仍是全长 12246
d_struct = 0.5*d_geo + 0.5*d_comp             # (12246, 12246)
logk_all = np.log10(meta["clarke"].values)    # 长度 12200
idx = rng.choice(n, size=4000)                # n = 12200
sub = d_struct[np.ix_(idx,idx)]               # 用过滤后下标索引未过滤矩阵
```
`figures/struct_vs_kL_dist.png` 上标的 Spearman 无效。

**bug C：merge 后长度未断言。**
`meta = meta.merge(elec, on="material_id", how="left")` 之后直接用 `soap_geo` 和 `meta` 按行号对齐。当前 `elec` 的 material_id 恰好唯一所以没炸，但这是运气 —— 一旦上游出现重复键，行会静默错位。必须加 `assert len(meta) == len(soap_geo)`。

## 问题 7 —— 统计检验顶在下限，等于没做【中等】

- 近邻重叠：`N_PERM=300` → 最小可报 p = 1/301 = 0.0033。表里**所有** p 都是 0.0033。
- Mantel：`N_PERM_MANTEL=100` → 最小 p = 1/101 = 0.0099。表里**所有** Mantel-p 都是 0.0099。

这两列不携带任何信息，只是在说"打到置换下限了"。
另外，Spearman 在 100 万个抽样对上做，n=12246 时任何 ρ 都会"显著" —— 显著性在这里没有意义，**效应量**才有。

**改法**：报置换 z 分数和效应量的 bootstrap 置信区间，不报 p。或者把置换次数提到能分辨的量级（≥10000）并明确报告"p < 1e-4"。

## 问题 8 —— 结论表内部自相矛盾【中等】

`view_distance_corr.csv` 里：

- `Elastic vs Structure` → Spearman **0.111**
- `Structure vs kL_clarke` → Spearman **0.362**
- `Elastic vs kL_clarke` → Spearman **0.650**

但 kL_clarke 是 Elastic 的确定性函数（R²=0.98）。如果 Structure 和 Elastic 只有 0.111 的关系，Structure 怎么可能和 Elastic 的函数有 0.362 的关系？

这个矛盾的来源就是问题 3 说的维度不对等 —— 高维↔高维（Elastic vs Structure）和高维↔1维（Structure vs kL）不是同一个统计量。**这条矛盾本身就是"当前度量不可比"的直接证据**，应该写进报告，而不是留在 CSV 里没人看。

## 问题 9 —— 混杂因素完全没有控制【中等】

Structure 视图 = 几何 SOAP + 组成 Hellinger。但组成**唯一决定**平均原子质量，几何+组成**唯一决定**密度。而 clarke 又是质量和密度的函数。

所以"Structure ↔ kL 有 8.2x 富集"至少一部分只是"结构信息里含有密度和质量信息"。

**你真正该问的是**：控制住密度、平均原子质量、弹性模量之后，结构还携带多少关于 κ_L 的**增量**信息？
这需要偏相关 / 条件互信息 / 嵌套模型消融，跨视图 kNN 重叠回答不了。

## 问题 10 —— SOAP 设置的两个隐患【中等】

- `species=["X"]` 抛弃全部化学信息，只留几何。这是刻意设计（README 有说），但要注意此时 SOAP 无法区分 NaCl 和 MgO —— 组成视图在补这个，0.5/0.5 权重却是**拍脑袋定的**，没有做敏感性分析。
- `average="inner"` 对所有格点做平均池化，把局部环境差异抹平了。
- **周期性 SOAP 对晶胞选择敏感**。MP 返回的 structure 有的是原胞有的是常规胞，你没有做标准化（`SpacegroupAnalyzer.get_primitive_standard_structure()`）。同一个材料换个胞的表示，SOAP 向量会变 —— 这是可复现性隐患。
- `d_geo / d_geo.max()` 的归一化被**单个最远离群点**控制。用 95 分位数归一化更稳。

## 问题 11 —— 工程卫生【轻度但要补】

缺：`requirements.txt`（带版本锁）、`LICENSE`、任何测试、数据版本记录（MP 数据库版本号、JARVIS 快照日期、下载日期）、随机种子的统一管理。

## 问题 12 —— 项目回答不了你想问的问题【方向性】

zT = S²σT / (κ_e + κ_L)。

当前项目只研究了 κ_L 的一个非晶极限近似，从未涉及 Seebeck 系数、电导率、功率因子、载流子浓度、温度。
即使上面 11 条问题全部修好、结论全部成立，它**仍然筛不出好的热电材料** —— 因为好的热电材料是 PF 高**且** κ_L 低的帕累托最优点，而你只看了一个轴的代理量。

`scripts/kl_verify_05_dual_channel_intersection.py` 是唯一朝这个方向走的（PF vs κ_L 帕累托），但只有 85 个样本，且文档里自己承认温度和数据来源不一致所以没算 zT。这个思路是对的，规模和严谨度不够。

---

# 第 2 部分：重构后的科学问题

把原来那个"结构相似是否 ⇒ κ_L 相似"改成三个可证伪、可回答、且直接服务于选材的问题：

**Q1（解耦性）**：控制住弹性/密度/质量这些已知的声子输运描述符后，晶体结构是否还携带关于真实 κ_L 的**增量**信息？
→ 方法：嵌套模型消融 + 偏相关 + 条件互信息。有标定曲线做参照。

**Q2（双通道正交性）**：电子输运品质（PF、S、σ）与晶格输运（κ_L）在描述符空间中是否可以独立调控？
→ 这才是热电材料设计的核心命题（"phonon-glass electron-crystal"）。
→ 方法：在同一批材料上同时拿到 PF 和 κ_L，看两者对同一组描述符的依赖是否正交。

**Q3（筛选）**：在 (PF, κ_L) 平面上的帕累托前沿是哪些材料？前沿材料在描述符空间里有什么共同特征？
→ 这直接产出候选清单，也是你说的"通过性质和参数之间的联系去探究优秀热电材料"。

---

# 第 3 部分：逐步提示词

---

## Step 0：环境与工程基线

**目标**：让仓库能被别人 clone 下来直接跑。

**验收**：`pip install -r requirements.txt && pytest` 通过；`grep -r "UxK6" .` 无结果（除 git 历史）。

```
我要重构一个材料科学数据分析仓库。先做工程基线，不要碰任何科学逻辑。

任务：
1. 删除 mp_kappaL/download_mp.py 里硬编码的 API key，改为：
   KEY = os.environ.get("MP_API_KEY")
   if not KEY: raise SystemExit("请设置环境变量 MP_API_KEY")
2. 修复 crossview_analysis.py 和 descriptor_analysis.py 里的 sys.path bug：
   它们 insert 的是 root/"jarvis_2d_te_atlas"/"scripts"，但 graph_utils.py 实际在 root/"scripts"。
   改成从 scripts/ 导入，或者把 graph_utils.py 做成包内模块。
3. 新建 requirements.txt，锁定版本（numpy pandas scipy scikit-learn dscribe pymatgen ase pyarrow matplotlib）。
   跑 pip freeze 拿实际版本号，不要凭记忆写。
4. 新建 LICENSE（MIT）。
5. 新建 config.py，集中管理所有魔法数字：K=10、随机种子、各种过滤阈值、SOAP 参数。
   所有脚本从这里读，不再各自硬编码。
6. 新建 tests/test_smoke.py：至少测 graph_utils 里的 hellinger_distance 和 soap_distance
   （对称性、非负性、对角为零、手算小例子）。

改完后把 requirements.txt 的实际内容和 pytest 输出给我看。
```

---

## Step 1：把已知的循环论证写成回归测试

**目标**：把"Clarke 是弹性量的函数"这件事从"隐藏的结论污染"变成"显式的正确性断言"。

**验收**：`tests/test_circularity.py` 通过，且 README 头条表里不再出现 Debye↔kL 的相关系数。

```
在 tests/ 下新建 test_circularity.py。

背景：Materials Project 的 thermal_conductivity.clarke 是 Clarke 最小热导率的解析公式，
κ_min = 0.87·k_B·M_a^(-2/3)·E^(1/2)·ρ^(1/6)，等价于 ∝ n^(2/3)·v_m。
它是 Debye 温度、声速、密度的确定性函数，不是独立的物理测量。

写测试：
1. 从 mp_kappaL/processed/views_meta.parquet 读数据
2. 剔除物理不可能值（见 Step 2 的清洗规则）
3. 断言 OLS 回归 log10(clarke) ~ log10(debye) + log10(density) 的 R² > 0.95
4. 断言 Spearman(debye, log10 clarke) > 0.95
5. 测试的 docstring 里明确写：这个测试通过 = 数据管线正确，
   不是物理发现；任何把 Debye↔clarke 的高相关当作结论的表述都是错的。

然后修改 README.md 和 reports/02_mp_kappaL_12246.md：
- 把 "Debye +0.99 / 纵声速 +0.96 / 剪切模量 +0.85" 这些从「结果」章节移到「流水线正确性检查」章节
- 把 "弹性(力学) ↔ kL 强联系（物理机制：声速/Debye）" 这一行改写为
  「弹性↔kL：数学恒等，用作流水线灵敏度校准，不构成物理结论」
- 在 README 顶部加一段显式说明：clarke/cahill 是非晶极限最小热导率，不是晶态 κ_L
```

---

## Step 2：源头数据清洗，可审计

**目标**：一次清洗，全局生效，被剔除的样本可查。

**验收**：`processed/rejected_samples.csv` 存在且非空；清洗后所有物理量落在合理区间。

```
重写 mp_kappaL/build_views.py 的清洗逻辑。

当前问题：views_meta.parquet 里存在物理不可能值，且过滤只在下游各脚本零散地做。
实测异常：debye 最大 3.48e8 K（46 条 >5000K）；shear_vrh 最大 5.6e13 GPa（57 条 >1000）；
v_long 最大 3.29e9 m/s（超光速，52 条 >2e4）；bulk_vrh 最大 1.37e8 GPa。

要求：
1. 在 config.py 里定义显式的物理合理性区间，每条都写注释说明依据：
   bulk_vrh:  (0.1, 700) GPa      # 金刚石 ~443 GPa，留余量
   shear_vrh: (0.1, 700) GPa
   debye:     (10, 3000) K        # 金刚石 ~2230 K
   v_long:    (500, 25000) m/s    # 金刚石纵波 ~18000 m/s
   v_trans:   (200, 20000) m/s
   density:   (0.3, 25) g/cm3
   clarke:    (0.01, 100) W/mK
   （这些是我给的起点，你要查证并在注释里写明每个上界的物理参照物）
2. 清洗在 build_views.py 里做一次，写成一个 clean_records(df) 函数
3. 输出 processed/rejected_samples.csv，列出每个被剔除的 material_id、
   触发的规则、以及该字段的原始值
4. 输出清洗前后的样本数和每条规则命中数到 stdout
5. 下游脚本（crossview_analysis, descriptor_analysis）删除所有自己的过滤逻辑，
   只读清洗后的数据
6. 加断言：清洗后所有字段都在区间内，否则 raise

跑完后把 rejected_samples.csv 的前 20 行和规则命中统计给我。
```

---

## Step 3：修复行对齐 bug 并加断言

**目标**：消灭静默的行错位。

**验收**：`tests/test_alignment.py` 通过。

```
修复行对齐 bug。

bug B（确认存在）：descriptor_analysis.py 里，meta 被 bulk_vrh<1000 过滤掉 46 行后，
仍然用全长的 comp_frac.npy 和 soap_geo.npy 计算 d_struct（12246×12246），
再用过滤后的下标（范围 0..12199）去索引它。figures/struct_vs_kL_dist.png 上的
Spearman 数值因此无效。

bug C（潜在）：crossview_analysis.py 里 meta.merge(elec, on="material_id", how="left")
之后直接按行号对齐 soap_geo，没有任何长度断言。

要求：
1. 建立统一约定：所有特征矩阵（soap_geo.npy、comp_frac.npy）必须与一个
   显式的 material_id 索引数组一起保存和加载。新增 processed/row_index.npy。
2. 写一个 load_aligned(names) 辅助函数：读入若干特征矩阵和 meta，
   按 material_id 做 inner join 并重排，返回严格对齐的结果。
   函数内部对每一步做 assert。
3. 所有脚本改用 load_aligned，禁止再出现裸的按行号对齐。
4. 新建 tests/test_alignment.py：构造一个故意错位的小例子，
   断言 load_aligned 能检测出来并抛异常。
5. 重跑 descriptor_analysis.py，对比修复前后 struct_vs_kL 的 Spearman 值，
   把两个数都报给我。

跑完后告诉我：修复前 Spearman 是多少，修复后是多少，差多少。
```

---

## Step 4：修复 kNN 的并列退化

**目标**：让 kNN 结果不依赖于输入行序。

**验收**：把输入行随机打乱后重跑，所有重叠指标变化 < 1%。

```
修复 kNN 并列退化。

实测：electronic_jarvis.parquet 里 gap_opt 有 69.65%（7687/11037）恰好等于 0
（金属的 OptB88 带隙）。这导致 Eg 视图的距离矩阵有一个 7687×7687 的全零块，
np.argpartition 的选择完全由数组存储顺序决定，没有任何物理含义。

注意：scripts/graph_utils.py 的 kNN_affinity 已经写了 tiebreak_seed 来处理这个问题，
但 crossview_analysis.py 用的是自己定义的 knn_sets 函数，没有用上。

要求：
1. 把 crossview_analysis.py 里的 knn_sets 删掉，统一用 graph_utils 里带
   tiebreak_seed 的实现
2. 更重要：加一个"并列诊断"步骤。对每个视图，报告：
   - 距离矩阵中恰好为零的非对角元占比
   - 每个点的第 k 近邻距离与第 (k+1) 近邻距离相等的比例
   把这个表存到 processed/knn_degeneracy.csv
3. 加一个稳健性测试：把输入行随机打乱 5 次，每次重跑全部跨视图重叠计算，
   报告各指标的变异系数。写到 processed/row_order_robustness.csv
4. 如果某个视图的零距离占比 > 20%，在报告里明确标注
   「该视图 kNN 结构由并列主导，重叠指标不可解释」

跑完后把 knn_degeneracy.csv 和 row_order_robustness.csv 给我。
```

---

## Step 5：重建电子性质视图（按 ID 而非化学式）

**目标**：消除跨库化学式匹配引入的伪重复。

**验收**：新的电子视图里，重复特征向量占比 < 1%（当前是 14.8%）。

```
重建电子性质视图。

当前问题（实测）：
- match_electronic.py 按 pymatgen reduced_formula 跨库匹配 MP 和 JARVIS，
  导致 14.8% 的 MP 材料共享同一个 canon：26 个 SiO2 多晶型、21 个 Al2O3、
  18 个 C、15 个 Si、15 个 SiC 各自被赋予完全相同的电子性质。
  这些人为重复点会直接抬高近邻重叠，制造假信号。
- 讽刺的是，download_mp.py 已经把 MP 自己的 band_gap 下到了
  raw/summary_bandgap.json（按 material_id 严格对应，零匹配噪声），但从未使用。

要求：
1. 废弃 match_electronic.py 的化学式匹配路径。保留文件但在顶部加
   DEPRECATED 说明和废弃原因。
2. 新写 mp_kappaL/build_electronic.py：
   - 从 raw/summary_bandgap.json 读 MP 自己的 band_gap、is_metal（按 material_id）
   - 加载 Ricci et al. BoltzTraP 数据库（见第 4 部分），它也是 MP-id 索引，
     提供 Seebeck、电导率/弛豫时间、功率因子、电子热导、电导有效质量
   - 严格按 material_id 做 join，绝不按化学式
   - 输出 processed/electronic_by_mpid.parquet
3. 加断言：material_id 唯一；重复特征向量（所有电子特征完全相同的行）占比 < 1%
4. 报告新旧两种匹配方式的样本数对比和重复率对比

跑完后告诉我：新方法拿到多少条电子数据，重复特征向量占比是多少。
```

---

## Step 6：接入真实 κ_L 数据（关键一步）

**目标**：把目标变量从"非晶极限公式"换成真实晶格热导率。

**验收**：新目标变量的 p1–p99 跨度 > 2 个数量级（当前 clarke 只跨 1.5 倍的四分位）。

```
接入真实晶格热导率数据。

背景：当前的目标变量 clarke/cahill 是非晶极限的最小热导率解析公式，
实测四分位区间只跨 2.4 倍（0.418–0.994 W/mK），且 log10(clarke) 可被
log(debye)+log(density) 以 R²=0.979 重构。它不是晶格热导率。

要求：搜索并接入至少两个独立的真实 κ_L 数据源，做成 processed/kappa_L_targets.parquet。
候选（你需要自己联网核实可获取性、许可、当前规模，不要相信我给的数字）：

  a) 实验：starrydata2（仓库里已有本地副本，N=137 交集）
     —— 权威但小，用作最终验证集
  b) 第一性原理 BTE：Togo 的 PhononDB / phono3py 数据集
     （103 个二元化合物那批），以及后续扩展的高通量 κ_L 数据集
  c) 近期的大规模自动化 κ_L 数据库（搜 "Phonix database auto-kappa lattice
     thermal conductivity" 和 "high-throughput lattice dynamics npj Comput Mater 2024"）
  d) AFLOW 的 AGL 模块（Slack 模型近似的 κ_L，规模较大但精度低于 BTE）
  e) 半 Heusler 专用数据集（~143 个），适合做单一结构族的对照

要求：
1. 每个源单独存一份，标注 method 列（experimental / BTE-phono3py / Slack-AGL / ...）
   和 temperature 列。绝不把不同方法的 κ_L 混在一起当同一个变量。
2. 尽量映射到 material_id；映射不上的保留结构文件路径
3. 输出数据源交叉表：每个源的 N、κ_L 的 p1/p50/p99、与其他源的重叠样本数
4. 对重叠样本，画不同源之间的 κ_L 一致性散点图（这本身就是有价值的结果）
5. 保留 clarke/cahill 作为「非晶下界参照」列，但明确不再作为主目标

跑完后把数据源交叉表和一致性散点图给我。如果找不到 N > 500 的真实 κ_L 数据，
停下来告诉我，我们改用「小样本高质量 + 大样本代理」的双轨设计。
```

---

## Step 7：建立度量标定曲线（这一步决定项目可不可信）

**目标**：知道跨视图 Spearman 读数 0.36 到底意味着多强的依赖。

**验收**：产出一条标定曲线图 + 一个 `calibrate(rho_observed) -> R²_estimated` 函数。

```
建立跨视图度量的标定曲线。这是整个重构里最重要的一步。

问题（实测）：Elastic 视图和 kL_clarke 在数学上是 R²=0.98 的确定性关系，
但你的跨视图 Spearman 只给出 +0.65。既然确定性关系的读数上限是 0.65，
那么 Structure 的 0.36 和 Eg 的 0.09 就无法解释 —— 我们不知道这把尺子的刻度。

另有一个内部矛盾佐证这一点：view_distance_corr.csv 里
Elastic vs Structure = 0.111，但 Structure vs kL_clarke = 0.362，
而 kL_clarke 是 Elastic 的确定性函数。这在逻辑上不自洽，
根源是高维↔高维和高维↔1维不是同一个统计量（测度集中效应）。

要求：
1. 写 mp_kappaL/calibration.py。用真实的 Structure 视图（不用合成特征，
   保留真实的距离分布结构），构造一系列合成目标变量：
   y_synth = f(真实描述符的某个线性/非线性组合) + noise
   调节噪声幅度，使 y_synth 与该组合的真实 R² 分别为
   0.0, 0.1, 0.2, ..., 0.9, 1.0
2. 对每个 R² 水平，跑完全相同的跨视图流程（kNN 重叠 + 距离 Spearman），
   记录读数。每个水平重复 20 次取均值和 95% 区间。
3. 画标定曲线：x 轴 = 真实 R²，y 轴 = 度量读数。分别为
   「高维↔1维」和「高维↔高维」两种情况画两条曲线。
4. 输出反函数 calibrate(rho_observed, dim_case) -> (R²_lo, R²_hat, R²_hi)
5. 用它去翻译现有结果：Structure↔kL 的 0.362 对应多大的真实 R²？
   置信区间是多少？Eg↔kL 的 0.094 能不能与 0 区分？

跑完后把标定曲线图和翻译后的结果表给我。
如果标定曲线显示该度量在 R² < 0.3 区间不可分辨，
我们就必须换度量（下一步会做）。
```

---

## Step 8：换成能回答增量信息问题的方法

**目标**：从"两个视图像不像"升级到"这个视图有没有独立贡献"。

**验收**：产出块消融表，每个描述符块的增量 R²（带 CV 误差棒）。

```
把跨视图 kNN 重叠替换/补充为能回答「增量信息」的方法。

背景：当前方法回答的是「视图 A 的近邻和视图 B 的近邻重不重合」，
它无法回答真正的科学问题：控制住已知描述符后，结构还携带多少关于 κ_L 的新信息？

而且存在明显混杂：Structure 视图 = 几何 SOAP + 组成 Hellinger，
而组成唯一决定平均原子质量、几何+组成唯一决定密度，
密度和质量又直接进入 κ 的表达式。所以「Structure↔κ 有关系」至少一部分是同义反复。

要求实现三种方法，都用嵌套交叉验证：
1. 【块消融】把描述符分成互斥的块：
   Block-C 组成（元素分数、平均原子质量、电负性统计等）
   Block-G 几何（SOAP、密度、空间群、配位数统计、原子堆积分数）
   Block-E 弹性（B、G、Debye、声速、泊松比、各向异性）
   Block-X 电子（band gap、Seebeck、PF、有效质量）
   用 GBM/随机森林预测 log κ_L，逐块加入和逐块剔除，
   报告每个块的增量 R²（5-fold CV，报均值±标准差）。
   **必须用按化学体系分组的 CV（GroupKFold by chemical system），
   不能用随机 CV** —— 否则同一族材料泄漏到测试集，R² 会虚高。
2. 【偏相关】对每个候选描述符，计算它与 log κ_L 的偏 Spearman，
   控制变量为 {密度, 平均原子质量, Debye}。
   这直接回答「除了已知的 n^(2/3)·v_m 之外还有什么」。
3. 【条件互信息】I(Structure ; κ_L | Elastic)，用 kNN 估计量（Kraskov 或 sklearn）。
   同时报 I(Structure ; κ_L) 作为对照，两者之差就是被弹性解释掉的部分。

输出 processed/block_ablation.csv、partial_corr.csv、conditional_mi.csv
以及一张「每个描述符块的增量贡献」的条形图（带误差棒）。

跑完后把三张表给我。特别关注：Block-G 在已有 Block-E 的情况下增量 R² 是多少。
```

---

## Step 9：修统计检验

**目标**：报效应量和置信区间，不报顶到下限的 p 值。

```
修复统计检验。

问题（实测）：view_overlap.csv 里所有 15 行的 p 都是 0.0033（= 1/301，置换下限）；
view_distance_corr.csv 里所有 Mantel-p 都是 0.0099（= 1/101，置换下限）。
这两列不携带信息。另外 Spearman 在 100 万抽样对上算，n=12246 时任何 ρ 都会「显著」，
显著性在这里无意义，效应量才有意义。

要求：
1. 删掉两个结果表里的 p 列，改为报告：
   - 置换 z 分数（保留）
   - 效应量的 bootstrap 95% 置信区间（1000 次 bootstrap，按 material_id 重采样）
   - 经 Step 7 标定曲线翻译后的「等效真实 R²」及其区间
2. 近邻重叠必须同时报绝对值，不能只报富集倍数。
   当前 Elastic↔kL 的 overlap=0.0205 意思是「平均 10 个近邻里共享 0.2 个」，
   而 "24.9x 富集" 是相对于 null=0.0008 的比值。
   富集倍数在 null 接近 0 时会剧烈放大，单独报是误导。
   报告里两个数必须并排出现。
3. 加 k 的敏感性分析：k ∈ {5, 10, 20, 50, 100}，看结论是否稳定。
   写到 processed/k_sensitivity.csv
4. 加多重比较校正（Benjamini-Hochberg），因为你在测十几个视图对。

跑完后把 k_sensitivity.csv 和新的结果表给我。
```

---

## Step 10：SOAP 与结构表征的稳健性

```
修复结构表征的可复现性隐患。

问题：
- MP 返回的 structure 有的是原胞有的是常规胞，代码没有做标准化。
  周期性 SOAP 对晶胞选择敏感 —— 同一材料换个表示，SOAP 向量会变。
- average="inner" 的平均池化抹平了局部环境差异。
- d_geo/d_geo.max() 的归一化被单个最远离群点控制。
- Structure = 0.5*几何 + 0.5*组成 的权重是拍脑袋定的，没做敏感性分析。

要求：
1. 在 build_views.py 里，所有结构先过
   SpacegroupAnalyzer(struct).get_primitive_standard_structure()
   做标准化。记录标准化前后原子数变化的统计。
2. 加一个 sanity check：对随机抽取的 100 个材料，
   分别用原胞和常规胞算 SOAP，报告向量余弦相似度的分布。
   标准化后应该都接近 1。
3. 归一化改用 95 分位数而非 max，并报告两种归一化下结论的差异。
4. 权重敏感性：w ∈ {0, 0.25, 0.5, 0.75, 1}（w=几何权重），
   全部跑一遍，画出结论随 w 的变化。写到 processed/weight_sensitivity.csv
5. 【可选，但推荐】增加一个对照表征：用 matminer 的
   ElementProperty + SineCoulombMatrix 或 XRD 图谱作为独立的结构表征，
   看主要结论是否表征无关。表征依赖的结论是弱结论。

跑完后把晶胞标准化的 sanity check 结果和 weight_sensitivity.csv 给我。
```

---

## Step 11：双通道分析（这才是热电材料的核心问题）

**目标**：回答 Q2 —— 电子输运和晶格输运能否独立调控。

```
建立双通道分析。这是项目真正对热电材料有价值的部分。

背景：zT = S²σT / (κ_e + κ_L)。好的热电材料需要 PF = S²σ 高、同时 κ_L 低。
「phonon-glass electron-crystal」这个概念的可检验版本就是：
存在一组材料，其电子输运品质由某组描述符控制，晶格输运由另一组描述符控制，
两组描述符可以独立调节。

当前仓库只有 scripts/kl_verify_05_dual_channel_intersection.py 做了一点点
（85 个样本的 PF vs κ_L 帕累托），思路对但规模和严谨度不够。

要求：
1. 用 Step 5 的 Ricci BoltzTraP 数据 + Step 6 的真实 κ_L，构建统一表：
   material_id | PF_n | PF_p | S_n | S_p | sigma_over_tau | kappa_e | kappa_L | T | doping
   注意 Ricci 数据是在给定掺杂浓度和温度下的，必须显式记录 T 和 doping，
   不能像 05 脚本那样把不同温度的量混在一起。
2. 分别做两次 Step 8 的块消融：
   目标 A = log(PF)，目标 B = log(κ_L)
   得到两张「描述符块 → 增量 R²」表。
3. 关键分析：计算两个目标各自的描述符重要性向量之间的相关性/夹角。
   如果接近正交 → 支持双通道解耦，意味着可以独立优化。
   如果高度相关 → 说明存在共同瓶颈，需要找出是什么。
4. 做条件分析：在固定 κ_L 分位区间内（比如 κ_L 的最低 20%），
   PF 的分布和描述符依赖是什么样？反之亦然。
   这直接回答「低 κ_L 是否必然牺牲 PF」。
5. 画出两个通道的描述符重要性对比图（双向条形图）。

输出 processed/dual_channel_importance.csv 和对应图表。
跑完后告诉我：两个通道的描述符重要性向量夹角是多少度。
```

---

## Step 12：筛选与候选清单

**目标**：产出你真正想要的东西 —— 一份有理有据的候选材料清单。

```
做帕累托筛选，产出候选材料清单。

要求：
1. 在 (PF, κ_L) 平面上计算帕累托前沿（max PF, min κ_L）。
   注意：PF 和 κ_L 必须在同一温度下。如果数据源温度不一致，
   要么插值到统一温度，要么分温度段分别做前沿，绝不混用。
   （原 05 脚本因为温度不一致而没算 zT，这个判断是对的 —— 现在把它做对。）
2. 计算估计 zT = S²σT / (κ_e + κ_L)。
   因为 Ricci 数据是常弛豫时间近似（σ/τ），必须：
   - 明确标注 τ 的假设值
   - 做 τ ∈ {1e-15, 1e-14, 1e-13} s 的敏感性分析
   - 报告 zT 的排序在不同 τ 下是否稳定（排序稳定比绝对值可靠）
3. 剔除不现实的候选：用 MP 的 energy_above_hull 过滤热力学不稳定相
   （阈值写在 config.py 里并说明理由），标注含毒性/稀缺元素（Pb, Cd, Hg, Te, Re...）。
4. 用 Step 8/11 的描述符重要性，分析前沿材料的共同特征：
   它们在哪些描述符上区别于非前沿材料？做一个对比统计表。
5. 对照已知的高性能热电体系（Bi2Te3、PbTe、SnSe、Mg3Sb2、方钴矿、
   Half-Heusler、笼状化合物）验证：你的方法能不能把它们排进前列？
   **这是最关键的验收标准。如果排不进去，方法有问题，不要急着报告新候选。**
6. 输出 reports/candidates.md：前 50 个候选，每个附上
   预测 PF、κ_L、估计 zT、稳定性、元素毒性/丰度标注、以及「为什么它上榜」的
   描述符归因（用 SHAP 或 permutation importance）。

跑完后先把第 5 步的验证结果给我 —— 已知高性能材料排在什么位置。
这个通过了再看候选清单。
```

---

## Step 13：重写报告

```
重写 README.md 和 reports/。

硬性要求：
1. 每一条结论必须标注三样东西：
   (a) 效应量 + 置信区间（不是 p 值）
   (b) 经 Step 7 标定后的等效真实依赖强度
   (c) 一句话的「这条结论在什么条件下会不成立」
2. 单独开一章「已知局限」，至少覆盖：
   - CRTA 假设对 σ 和 PF 的影响
   - 不同 κ_L 数据源之间的系统偏差
   - DFT 带隙低估对电子输运的影响
   - 训练分布之外的外推风险
   - 描述符块之间无法完全解耦的部分
3. 单独开一章「本项目 v0 版本的错误及修正」，诚实列出：
   - clarke 是最小热导率而非 κ_L，导致 v0 的主结论是循环论证
   - 度量未标定，导致 v0 的强/中/弱判断无依据
   - 化学式跨库匹配导致 14.8% 的伪重复
   - 三个行对齐/导入 bug
   写清楚每条的修正方式和修正后结论的变化。
   **这一章会让项目可信度上升而不是下降。** 主动披露自己发现的错误
   是研究工作最有价值的部分之一。
4. 结论表里，任何一个数字如果没有对应的验收测试，就不能出现在表里。
5. 加「如何复现」章节：从空环境到出图的完整命令序列，
   标注每步的预期运行时间和内存需求。
6. 加数据版本记录：MP 数据库版本、JARVIS 快照日期、各数据集的下载日期和 DOI。

写完后把新 README 的目录结构给我看。
```

---

# 第 4 部分：数据源清单

执行 Step 5 和 Step 6 时需要，**每一个都要自己联网核实当前可获取性、许可和规模**，不要相信下面的数字。

## 电子输运

| 源 | 内容 | 索引 | 备注 |
|---|---|---|---|
| **Ricci et al. 2017, Sci Data** | 约 48000 个材料的 BoltzTraP 电子输运：S、σ/τ、κ_e、电导有效质量 | MP material_id | Dryad `doi:10.5061/dryad.gn001`；figshare 上另有整理好的表格版 `ricci_boltztrap_mp_tabular`，给出 10¹⁸ cm⁻³ / 300 K 下 n 型和 p 型的均值，以及在温度 [100,1300] K、掺杂 [10¹⁶,10²¹] cm⁻³ 范围内的最优值。**这个表格版最适合直接用。** 常弛豫时间近似，τ 依赖量要谨慎 |
| **JARVIS-DFT** | 约 36000 个体材料的 BoltzTraP 结果 | JVASP id | Choudhary et al. 的热电数据驱动工作，可与 MP 交叉验证 |
| **MP summary** | band_gap、is_metal、efermi | MP material_id | 你已经下载了但没用 |

## 晶格热导率

| 源 | 方法 | 规模量级 | 备注 |
|---|---|---|---|
| **starrydata2** | 实验 | 你已有交集 137 | 权威，作最终验证集 |
| **PhononDB / phono3py（Togo）** | BTE 有限位移 | 百量级二元化合物 | 高质量，一致的计算设置 |
| **Phonix / auto-kappa 类数据库** | Peierls-BTE + Wigner 相干 | 较大 | 近期工作，搜 "Phonix database lattice thermal conductivity"，含约 1900 条 κ_L < 1 W/mK 的低热导材料 |
| **AFLOW AGL** | Slack 模型近似 | 万量级 | 规模大但精度低于 BTE，只能作辅助 |
| **半 Heusler 专用集** | BTE | 百量级 | 单一结构族，适合做「控制结构族」的对照实验 |

**重要**：不同方法的 κ_L 不能直接混用。BTE 和 Slack 模型的系统偏差可能超过一个数量级。
Step 6 要求你先画一致性散点图，就是为了量化这个偏差。

## 稳定性与筛选

- MP `energy_above_hull`、`formation_energy_per_atom`
- 元素丰度 / 价格 / 毒性：pymatgen 的 `Element` 有部分属性，其余需自建查找表

---

# 第 5 部分：重构后 README 应该长什么样

```markdown
# 热电材料的双通道描述符分析

## 研究问题
Q1 控制已知声子输运描述符后，晶体结构对真实 κ_L 是否有增量信息？
Q2 电子输运品质与晶格输运在描述符空间中是否正交？
Q3 (PF, κ_L) 帕累托前沿上有哪些材料，它们的共同特征是什么？

## 一句话结论
[每条都带效应量、CI、标定后的真实依赖强度、以及失效条件]

## 方法可信度
- 正对照：[弹性→κ_min 的确定性关系，用于校准流水线灵敏度]
- 负对照：[随机重排标签，度量应回到基线]
- 标定曲线：[度量读数 → 真实依赖强度 的映射，见 figures/calibration.png]
- 已知高性能热电体系的召回率：[Bi2Te3/PbTe/SnSe/... 在候选排名中的位置]

## 数据
[每个源：规模、方法、温度、许可、下载日期、DOI]
[数据源之间的一致性分析]

## 结果
[Q1/Q2/Q3 各自的结果，全部带误差棒]

## 已知局限
[至少 5 条，具体到会影响哪个结论]

## v0 版本的错误与修正
[诚实列出，这一章提升可信度]

## 复现
[完整命令序列 + 预期耗时 + 内存需求]
```

---

# 附：优先级建议

如果时间有限，按这个顺序做，每一步都独立有价值：

1. **Step 0**（吊销 key + 修 import bug）—— 半小时，必须做
2. **Step 1**（把循环论证写成测试 + 改 README 表述）—— 一小时，这一步就能让项目从"结论错误"变成"诚实"
3. **Step 2、3、4**（清洗 + 对齐 + 并列）—— 一天，修掉所有已确认的 bug
4. **Step 7**（标定曲线）—— 这是让现有结论变得可解释的最小代价方案。做完这一步，你甚至可以不换数据源，只是把结论重新表述为标定后的区间
5. **Step 6 + 8**（真实 κ_L + 增量信息方法）—— 这两步做完，Q1 才真正被回答
6. **Step 11 + 12**（双通道 + 筛选）—— 这才是「探究优秀热电材料」

Step 1 单独做完就已经是一次有意义的提交。不要等到全部做完才推。
