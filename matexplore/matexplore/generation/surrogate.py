# -*- coding: utf-8 -*-
"""替代模型：加载训练好的 Ridge 系数 + 标准化器，对 14 维易算特征打分。"""
import json
import os
import numpy as np

from ..hypothesis.cheap_features import feature_vector, FEATURE_ORDER

MODEL_NAMES = ["n_ZT_e", "p_ZT_e", "n_log10_PF", "p_log10_PF"]


class SurrogateScorer:
    def __init__(self, models_dir):
        self.models = {}
        for name in MODEL_NAMES:
            p = os.path.join(models_dir, name + ".json")
            if os.path.exists(p):
                with open(p) as f:
                    self.models[name] = json.load(f)

    def predict(self, feat_dict):
        """返回 {n_ZT_e, p_ZT_e, n_log10_PF, p_log10_PF} 预测值。"""
        x = feature_vector(feat_dict)
        out = {}
        for name, blob in self.models.items():
            order = blob["features"]
            idx = [FEATURE_ORDER.index(f) for f in order]
            xm = x[idx].copy()
            xm = (xm - np.array(blob["feature_mean"])) / np.array(blob["feature_scale"])
            xm = np.nan_to_num(xm, nan=0.0)
            out[name] = float(np.dot(xm, np.array(blob["coef"])) + blob["intercept"])
        return out

    def merit(self, pred, carrier="p", weight_zt=0.7, weight_pf=0.3):
        """综合热电优值打分：ZT_e 与 PF 的加权(均为越大越好)。"""
        zt = pred.get(f"{carrier}_ZT_e", np.nan)
        pf = pred.get(f"{carrier}_log10_PF", np.nan)
        return weight_zt * zt + weight_pf * pf

