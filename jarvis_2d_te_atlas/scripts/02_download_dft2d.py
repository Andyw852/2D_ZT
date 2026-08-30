"""Phase B/C - Step 2: 下载 JARVIS dft_2d。

数据源说明（重要）：
- jarvis-tools 默认下载源 get_db_info()["dft_2d"][0] = https://ndownloader.figshare.com/files/38521268
- 在当前执行环境中，figshare.com / api.figshare.com / ndownloader.figshare.com 全部返回 HTTP 403
  （AWS ALB 网络层阻断，非数据文件问题），已在 reports/download_probe.json 记录。
- 因此改用官方 NIST 托管的 JARVIS-DFT OPTIMADE API 作为权威数据源：
    https://jarvis.nist.gov/optimade/jarvisdft/v1/structures/
  该 API 的 entry id 前缀为 "dft_2d_"（对应 dft_2d 数据集），_jarvis_source == "dft_2d"，
  且包含完整热电性质（n/p-Seebeck, n/p-powerfact, n/pcond, n/pkappa, 有效质量等）。
- OPTIMADE 端点为 2D + 3D 混合；通过 id 前缀 "dft_2d_" 精确筛选 2D 数据。
- 分页：page 参数从 1 递增，每页固定 20 条（page_limit 被该实现忽略）。

本脚本下载全部 dft_2d 记录并保存快照。
"""
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import requests

BASE = "https://jarvis.nist.gov/optimade/jarvisdft/v1/structures/"
HEADERS = {"User-Agent": "Mozilla/5.0 (JARVIS dft_2d research download)"}

root = Path(__file__).resolve().parents[1]
raw_dir = root / "data" / "raw" / "jarvis"
raw_dir.mkdir(parents=True, exist_ok=True)
reports = root / "reports"
reports.mkdir(exist_ok=True)

records = []
page = 1
consecutive_empty = 0
max_pages = 300
seen_3d = False
while page <= max_pages:
    r = requests.get(BASE, params={"filter": "nelements>0", "page": page, "page_limit": 20},
                     timeout=120, headers=HEADERS)
    r.raise_for_status()
    d = r.json()
    data = d.get("data", [])
    n_2d = 0
    n_3d = 0
    for entry in data:
        eid = entry.get("id", "")
        if eid.startswith("dft_2d_"):
            records.append(entry)
            n_2d += 1
        elif eid.startswith("dft_3d_"):
            n_3d += 1
    if n_3d > 0:
        seen_3d = True
    print(f"page={page}: returned={len(data)} dft_2d={n_2d} dft_3d={n_3d} cumulative_2d={len(records)}")
    if n_2d == 0:
        consecutive_empty += 1
    else:
        consecutive_empty = 0
    # 停止条件：已经看到 3D 数据（说明 2D 块结束），且至少连续 3 页没有新的 2D
    if seen_3d and consecutive_empty >= 3 and len(records) > 500:
        print("Stopping: passed the dft_2d block.")
        break
    # 没有 next 链接时停止
    nxt = d.get("links", {}).get("next")
    if not nxt:
        print("No next link; stopping.")
        break
    page += 1
    time.sleep(0.05)

print("\nTOTAL dft_2d records:", len(records))

snapshot_path = raw_dir / "dft_2d_snapshot.json"
with snapshot_path.open("w", encoding="utf-8") as f:
    json.dump(records, f, ensure_ascii=False)

meta = {
    "execution_date": datetime.now(timezone.utc).isoformat(),
    "source": "JARVIS-DFT OPTIMADE API (official NIST host)",
    "base_url": BASE,
    "default_figshare_url_blocked": "https://ndownloader.figshare.com/files/38521268 (HTTP 403 in this environment)",
    "record_count": len(records),
    "pages_fetched": page,
    "saved_to": str(snapshot_path),
}
(reports / "download_metadata.json").write_text(
    json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8"
)
print(json.dumps(meta, indent=2, ensure_ascii=False))
