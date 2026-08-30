"""Phase G: 从官方 NIST static XML server 恢复原始 3 个输运本征值。

数据源：https://www.ctcms.nist.gov/~knc6/static/JARVIS-DFT/{JID}（官方 JARVIS-DFT 静态文件）
每个 JID 返回完整 XML，其中 <boltztrap_info> 含 p/n 的 seeb/cond/pf/kappa 三个本征值，
<effective_mass> 含 electron/hole_mass_300K 三个本征值。
"""
import json
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests

BASE = "https://www.ctcms.nist.gov/~knc6/static/JARVIS-DFT/"
HEADERS = {"User-Agent": "Mozilla/5.0 (JARVIS transport eigenvalue recovery)"}

root = Path(__file__).resolve().parents[1]
snapshot = root / "data" / "raw" / "jarvis" / "dft_2d_snapshot.json"
recs = json.loads(snapshot.read_text(encoding="utf-8"))
jids = [r["attributes"]["_jarvis_jid"] for r in recs]
print("total JIDs:", len(jids))

BT_TAGS = ["pseeb", "pcond", "ppf", "pkappa", "nseeb", "ncond", "npf", "nkappa"]
MASS_TAGS = ["electron_mass_300K", "hole_mass_300K"]

def _unquote(s):
    return s.strip().strip("'").strip('"').strip()

def _parse_inner(inner):
    out = {}
    for tag in BT_TAGS + MASS_TAGS:
        m = re.search(r"<" + tag + r">(.*?)</" + tag + r">", inner, re.S)
        if m:
            v = _unquote(m.group(1))
            vals = [x.strip() for x in v.split(",")]
            out[tag] = vals
    return out

def parse_xml(text):
    out = {}
    m = re.search(r"<boltztrap_info>(.*?)</boltztrap_info>", text, re.S)
    if m:
        inner = _unquote(m.group(1))
        if inner:
            out.update(_parse_inner(inner))
    m2 = re.search(r"<effective_mass>(.*?)</effective_mass>", text, re.S)
    if m2:
        inner = _unquote(m2.group(1))
        if inner:
            out.update(_parse_inner(inner))
    return out

def fetch_one(jid):
    url = BASE + jid
    try:
        r = requests.get(url, timeout=90, headers=HEADERS)
        if r.status_code != 200:
            return jid, None, f"HTTP_{r.status_code}"
        eig = parse_xml(r.text)
        return jid, eig, "OK"
    except Exception as e:  # noqa
        return jid, None, f"ERR_{type(e).__name__}"

results = {}
errors = {}
t0 = time.time()
with ThreadPoolExecutor(max_workers=12) as ex:
    futs = {ex.submit(fetch_one, j): j for j in jids}
    done = 0
    for fut in as_completed(futs):
        jid, eig, status = fut.result()
        done += 1
        if status == "OK":
            results[jid] = eig
        else:
            errors[jid] = status
        if done % 100 == 0:
            print(f"  {done}/{len(jids)}  ok={len(results)}  err={len(errors)}  elapsed={time.time()-t0:.0f}s")

# 保存原始提取结果
out_json = root / "data" / "processed" / "transport_eigenvalues_raw.json"
out_json.parent.mkdir(parents=True, exist_ok=True)
out_json.write_text(json.dumps(results, ensure_ascii=False), encoding="utf-8")

# 统计
n_bt = sum(1 for v in results.values() if "nseeb" in v or "pseeb" in v)
n_mass = sum(1 for v in results.values() if "electron_mass_300K" in v)
print(f"\nDONE: recovered={len(results)} errors={len(errors)}")
print(f"with boltztrap (any seeb): {n_bt}")
print(f"with effective_mass: {n_mass}")
print(f"saved: {out_json}")
if errors:
    print("error counts:", {k: sum(1 for v in errors.values() if v.startswith(k)) for k in sorted(set(errors.values()))})
