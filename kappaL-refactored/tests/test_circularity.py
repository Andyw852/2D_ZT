"""Step 1：把已知的循环论证写成回归测试。

Clarke 最小热导率的定义式是 κ_min = 0.87·k_B·M_a^(-2/3)·E^(1/2)·ρ^(1/6)，
等价于 κ ∝ n^(2/3)·v_m；Debye 温度 Θ_D ∝ n^(1/3)·v_m。
二者共享同一组变量（声速 v_m、数密度 n、质量、密度），所以
log10(clarke) 能被 log(debye)+log(density) 以 R²≈0.98 重构，Spearman(debye, log κ)≈0.99。

**这个测试通过 = 数据管线正确（能恢复已知的解析关系），不是物理发现。**
任何把 Debye↔clarke 的高相关当作「物理结论」的表述都是错的 —— 它是代数恒等式。
"""
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import config
from mp_kappaL.data_utils import clean_records


def _load_clean():
    meta = pd.read_parquet(config.PROC_DIR / "views_meta.parquet")
    clean, _ = clean_records(meta)
    return clean.dropna(subset=["clarke", "debye", "density"]).reset_index(drop=True)


def test_clarke_is_algebraic_function_of_debye_and_density():
    d = _load_clean()
    y = np.log10(d["clarke"].values.astype(float))
    X = np.column_stack([np.ones(len(d)),
                         np.log(d["debye"].values.astype(float)),
                         np.log(d["density"].values.astype(float))])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    r2 = 1.0 - resid.var() / y.var()
    assert r2 > 0.95, f"log10(clarke) ~ log(debye)+log(density) 的 R²={r2:.4f} 未达到 0.95"


def test_spearman_debye_clarke_near_one():
    d = _load_clean()
    y = np.log10(d["clarke"].values.astype(float))
    rho, _ = spearmanr(d["debye"].values.astype(float), y)
    assert rho > 0.95, f"Spearman(debye, log10 clarke)={rho:.4f} 未达到 0.95"


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn()
        print(f"PASS {fn.__name__}")
    print(f"all {len(fns)} tests passed")
