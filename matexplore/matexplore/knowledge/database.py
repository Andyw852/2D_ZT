# -*- coding: utf-8 -*-
"""知识库读取：既有数据库 + 派生特征面板 + 相关面板 + 替代模型。

对外暴露:
  load_panel(cfg)        -> DataFrame(1103 行, jid + 特征 + n/p 目标)
  load_correlation(cfg)  -> DataFrame(特征 x 目标 x carrier 的 Spearman)
  load_surrogate(cfg)    -> dict { "{carrier}_{target}": model_blob }
"""
import json
import os
import pandas as pd


def load_panel(cfg):
    return pd.read_csv(cfg.knowledge.panel_csv)


def load_correlation(cfg):
    return pd.read_csv(cfg.knowledge.correlation_csv)


def load_surrogate(cfg):
    d = cfg.knowledge.surrogate_models_dir
    models = {}
    for fn in os.listdir(d):
        if fn.endswith(".json"):
            with open(os.path.join(d, fn)) as f:
                models[fn[:-5]] = json.load(f)
    return models


def load_snapshot(cfg):
    """加载 JARVIS dft_2d 快照(列表 of {id, type, attributes})。"""
    with open(cfg.knowledge.structure_source) as f:
        return json.load(f)
