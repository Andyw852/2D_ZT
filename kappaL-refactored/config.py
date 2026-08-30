"""集中管理所有魔法数字、路径与物理合理性区间。

所有脚本从这里读配置，不再各自硬编码（Step 0 / Step 2 / Step 7 / Step 12 要求）。
"""
from pathlib import Path

# ---- 路径 ----
ROOT = Path(__file__).resolve().parent          # kappaL-refactored/
MP_DIR = ROOT / "mp_kappaL"
RAW_DIR = MP_DIR / "raw"
PROC_DIR = MP_DIR / "processed"
FIG_DIR = MP_DIR / "figures"
SCRIPTS_DIR = ROOT / "scripts"
REPORTS_DIR = ROOT / "reports"
TESTS_DIR = ROOT / "tests"

# 外部大数据集（starrydata2、JARVIS dft_3d）体积大且可再生，存放在工作区而非本仓库内。
# 路径：<workspace>/jarvis_2d_te_atlas/data/raw/external/
EXTERNAL_DATA_DIR = ROOT.parent / "jarvis_2d_te_atlas" / "data" / "raw" / "external"

# ---- 随机种子（统一管理，Step 0 要求）----
SEED = 42
N_PERM = 200               # 近邻重叠置换次数（估计 z 的 null std；Step 9 报 z+CI 不需顶到下限的 p）
N_PERM_MANTEL = 50         # Mantel 置换次数（同上）
N_BOOTSTRAP = 500          # 效应量 bootstrap 次数（Step 9）
N_ORDER_SHUFFLES = 5       # 行序稳健性重排次数（Step 4）
ROW_ORDER_AUDIT_N = 3000   # 行序是离散不变性检查；固定子样本避免重复复制 12k² 距离矩阵
WEIGHT_AUDIT_N = 3000      # SOAP 权重/归一化敏感性固定子样本
N_CALIB_REPEAT = 20        # 标定曲线每个 R² 水平的重复次数（Step 7）
N_CALIB_SAMPLE = 1500      # 标定曲线每次抽样点数（控制计算量，Step 7）

# ---- 跨视图度量 ----
K = 10                     # kNN 近邻数
M_SAMP = 1_000_000         # 成对距离 Spearman 抽样上限
K_SENSITIVITY = [5, 10, 20, 50, 100]   # k 敏感性分析（Step 9）
BH_ALPHA = 0.05            # Benjamini-Hochberg 多重比较 FDR（Step 9）

# ---- SOAP 参数（Step 10）----
SOAP_SPECIES = ["X"]       # dummy 物种 = 几何-only（有意为之，组成由 Hellinger 单独表征）
SOAP_R_CUT = 6.0
SOAP_N_MAX = 6
SOAP_L_MAX = 6
SOAP_SIGMA = 1.0
SOAP_AVERAGE = "inner"

# ---- 结构视图权重（Step 10 做敏感性扫描 w ∈ {0,0.25,0.5,0.75,1}）----
W_GEO = 0.5                # 几何 SOAP 权重
W_COMP = 0.5               # 组成 Hellinger 权重
W_GEO_GRID = [0.0, 0.25, 0.5, 0.75, 1.0]

# ---- 物理合理性区间（Step 2；每条注明物理参照物）----
# 起点来自审计文档实测（debye 最大 3.48e8 K、v_long 超光速、bulk_vrh 1.37e8 GPa），
# 上界以金刚石为最硬参照物并留余量。
PHYSICAL_RANGES = {
    "bulk_vrh":  (0.1, 700.0),     # GPa；金刚石体模量 ~443 GPa，留余量到 700
    "shear_vrh": (0.1, 700.0),     # GPa；金刚石剪切模量 ~535 GPa，留余量
    "debye":     (10.0, 3000.0),   # K；金刚石 Debye ~2230 K，留余量到 3000
    "v_long":    (500.0, 25000.0), # m/s；金刚石纵波 ~18000 m/s，留余量（真空光速 3e8 为绝对上限）
    "v_trans":   (200.0, 20000.0), # m/s；金刚石横波 ~12000 m/s，留余量
    "density":   (0.3, 25.0),      # g/cm3；锇 ~22.6 g/cm3，留余量
    "clarke":    (0.01, 100.0),    # W/mK；非晶极限 κ_min 的合理取值域
    "snyder_acoustic": (0.001, 5000.0), # W/mK；300 K Snyder 解析估计，覆盖金刚石量级并剔除数值爆炸
}

# ---- 稳定性 / 筛选（Step 12）----
ENERGY_ABOVE_HULL_MAX = 0.05       # eV/atom；热力学稳定/近稳定相阈值
TAU_VALUES = [1e-15, 1e-14, 1e-13] # s；CRTA 弛豫时间敏感性分析（zT 排序稳定性）
TOXIC_SCARCE_ELEMENTS = ["Pb", "Cd", "Hg", "Te", "Re", "Tl", "As", "Bi", "Sb"]

# ---- 电子视图（Step 5）----
MAX_DUP_FEATURE_FRACTION = 0.01    # 重复特征向量占比阈值（验收 < 1%）

# ---- 视图定义（Step 8 块消融的互斥描述符块）----
BLOCKS = {
    "C": ["composition", "avg_mass", "num_elements"],                         # 组成块
    "G": ["soap"],                                                              # 几何-only 块
    "E": ["bulk_vrh", "shear_vrh", "debye", "v_long", "v_trans", "density", "nsites"],
    "X": ["band_gap", "is_metal", "efermi"],                                 # 电子块
}
# 偏相关 / 条件互信息里的已知混杂控制变量（Step 8）
CONFOUNDERS = ["density", "avg_mass", "debye"]
