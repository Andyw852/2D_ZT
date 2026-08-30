#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""用 MatGL 预训练模型预测带隙。在 3090 matgl_venv 里跑；经 hf-mirror 下载权重。"""
import os, sys, json, warnings
warnings.filterwarnings("ignore")
os.environ["HF_ENDPOINT"] = "https://hf-mirror.com"

from pymatgen.core import Structure
import matgl

MODEL = "MEGNet-BandGap-mfi-MP-2019.4.1"

def main(paths):
    model = matgl.load_model(MODEL)
    out = {}
    for p in paths:
        try:
            s = Structure.from_file(p)
            g = model.predict_structure(s)
            v = float(g) if not hasattr(g, "__len__") else float(g[0])
            out[p] = round(v, 5)
        except Exception as e:
            out[p] = {"error": str(e)[:200]}
    print(json.dumps(out, ensure_ascii=False))

if __name__ == "__main__":
    main(sys.argv[1:])
