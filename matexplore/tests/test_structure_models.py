# -*- coding: utf-8 -*-
"""回归测试：make_poscar 的笛卡尔→分数坐标转换（JARVIS 行存储晶格）。

历史 bug：frac = coords @ inv(lat.T) 把 2H 结构分数坐标从 (1/3,2/3) 错成
(0.7887,-0.122)，W-Se 键长 2.552→2.369 Å。正确应为 inv(lat)。
"""
import sys
import unittest
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "matexplore" / "generation"))
from structure_models import make_poscar  # noqa: E402

# WSe2 2H（JARVIS dft_2d，JVASP）真实 lattice_vectors + cartesian_site_positions
LATTICE = np.array([
    [1.662139, -2.878908, 0.0],
    [1.662139, 2.878908, 0.0],
    [0.0, 0.0, 35.068951],
], dtype=float)
SPECIES = ["Se", "Se", "W"]
COORDS = np.array([
    [1.66214, 0.959639, 2.085145],
    [1.66214, 0.959639, 5.449314],
    [1.66214, -0.959639, 3.767247],
], dtype=float)


def _parse_direct(poscar_text):
    lines = poscar_text.splitlines()
    i = next(i for i, l in enumerate(lines) if l.strip() == "Direct")
    return np.array([[float(t) for t in l.split()[:3]] for l in lines[i + 1:]])


class TestMakePoscar(unittest.TestCase):
    def test_wse2_fractional_coords(self):
        """2H WSe2 的 W 应在 (2/3,1/3)、Se 应在 (1/3,2/3)。"""
        frac = _parse_direct(make_poscar(LATTICE, SPECIES, COORDS))
        self.assertTrue(np.allclose(frac[2, :2], [2 / 3, 1 / 3], atol=1e-3), frac)
        self.assertTrue(np.allclose(frac[0, :2], [1 / 3, 2 / 3], atol=1e-3), frac)

    def test_cartesian_roundtrip(self):
        """Direct 坐标乘晶格应还原原始笛卡尔坐标（行存储）。"""
        frac = _parse_direct(make_poscar(LATTICE, SPECIES, COORDS))
        cart = frac @ LATTICE
        self.assertTrue(np.allclose(cart, COORDS, atol=1e-3))

    def test_not_transposed(self):
        """显式对比：inv(lat) 与 inv(lat.T) 给出不同结果，且只有 inv(lat) 正确。"""
        lat = LATTICE
        coords = COORDS
        frac_ok = coords @ np.linalg.inv(lat)
        frac_bad = coords @ np.linalg.inv(lat.T)
        self.assertFalse(np.allclose(frac_ok, frac_bad, atol=1e-6))
        self.assertTrue(np.allclose(frac_ok @ lat, coords, atol=1e-3))


if __name__ == "__main__":
    unittest.main()
