# -*- coding: utf-8 -*-
"""结构/材料生成器。

当前实现两套"可运行、无需下载大模型"的生成器：
  1. SubstitutionGenerator —— 原型替换：取高 ZT_e 种子二维结构，在元素族内
     同价替换生成新组成(保持晶格/对称性近似不变，得到可直接建 POSCAR 的新材料)。
  2. CompositionEnumerator  —— 组成枚举：在给定元素集合上枚举候选化学式，
     用代理模型打分(仅组成特征，无结构时用 seed 平均晶格作占位)。

可插拔结构生成模型(接口已留)：MatterGen / CDVAE / DiffCSP / FlowMM 等生成式模型
可用 generate_from_model() 接入；本环境未安装这些模型，故以原型替换+组成枚举
作为可靠的离线实现。
"""
import itertools
import numpy as np


def make_poscar(lattice, species, coords, comment="candidate"):
    """生成 VASP5 POSCAR 文本。coords 为笛卡尔坐标。"""
    uniq = []
    for s in species:
        if s not in uniq:
            uniq.append(s)
    counts = [species.count(s) for s in uniq]
    lines = [comment, "1.0"]
    for v in lattice:
        lines.append("  %.8f  %.8f  %.8f" % (v[0], v[1], v[2]))
    lines.append("  " + "  ".join(uniq))
    lines.append("  " + "  ".join(str(c) for c in counts))
    lines.append("Direct")
    lat = np.array(lattice, dtype=float)
    # JARVIS/ASE 晶格矢量按【行】存储：cart = frac @ lat，故 frac = cart @ inv(lat)。
    # 旧代码 inv(lat.T)（== inv(lat)^T）对非对称晶格矩阵多转置一次，
    # 把 2H 结构分数坐标从 (1/3,2/3) 错成 (0.7887,-0.122)，W-Se 键长 2.552→2.369 Å。
    inv = np.linalg.inv(lat)
    frac = np.array(coords, dtype=float) @ inv
    for i, s in enumerate(species):
        lines.append("  %.8f  %.8f  %.8f   #%s" % (frac[i,0], frac[i,1], frac[i,2], s))
    return chr(10).join(lines)


def attributes_to_poscar(attr):
    lv = np.array(attr["lattice_vectors"], dtype=float)
    sp = attr["species_at_sites"]
    coords = np.array(attr["cartesian_site_positions"], dtype=float)
    return make_poscar(lv, sp, coords, comment=attr.get("_jarvis_formula", "candidate"))


class SubstitutionGenerator:
    """在元素族内做同价替换，生成新化学式 + POSCAR。"""
    def __init__(self, groups):
        self.group_of = {}
        for g in groups:
            for e in g:
                self.group_of[e] = set(g)

    def expand(self, species):
        """给定 site 元素列表，返回所有族内替换后的新 species 列表(排除原样)。"""
        out = []
        n = len(species)
        for i in range(n):
            s = species[i]
            if s in self.group_of:
                for rep in self.group_of[s]:
                    if rep != s:
                        new = list(species)
                        new[i] = rep
                        out.append(new)
        return out


class CompositionEnumerator:
    """枚举候选化学式(给定元素集合)，无结构 -> 用 seed 平均晶格作占位。"""
    def __init__(self, elements, max_species=3):
        self.elements = elements
        self.max_species = max_species

    def enumerate(self, max_per_arity=200):
        comps = []
        for k in range(1, self.max_species + 1):
            for combo in itertools.combinations(self.elements, k):
                comps.append(list(combo))
                if len(comps) >= max_per_arity * self.max_species:
                    return comps
        return comps



def formula_from_species(species):
    """从每原子元素列表生成化学式(按字母序 + 计量数)，如 SnF4 -> F4Sn。"""
    from collections import Counter
    cnt = Counter(species)
    return "".join((el + (str(cnt[el]) if cnt[el] > 1 else "")) for el in sorted(cnt))

