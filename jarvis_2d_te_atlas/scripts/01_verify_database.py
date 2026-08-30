"""Phase B - Step 1: 验证 JARVIS dft_2d 真实官方数据源。

禁止硬编码 FigShare URL；从当前安装的 jarvis-tools 实时读取官方配置。
"""
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import jarvis

from jarvis.db.figshare import get_db_info

root = Path(__file__).resolve().parents[1]
reports = root / "reports"
reports.mkdir(exist_ok=True)

info = get_db_info()

if "dft_2d" not in info:
    raise RuntimeError("Current jarvis-tools does not provide dft_2d.")

dataset_info = info["dft_2d"]

print("Download URL:", dataset_info[0])
print("Internal JSON:", dataset_info[1])
print("Description:", dataset_info[2])
print("Reference:", dataset_info[3])

record = {
    "execution_date": datetime.now(timezone.utc).isoformat(),
    "python_version": sys.version,
    "jarvis_tools_version": getattr(jarvis, "__version__", "unknown"),
    "dataset": "dft_2d",
    "download_url": dataset_info[0],
    "internal_json": dataset_info[1],
    "description": dataset_info[2],
    "reference": dataset_info[3],
    "all_datasets_available": sorted(info.keys()),
}

out = reports / "database_source.txt"
with out.open("w", encoding="utf-8") as f:
    for k, v in record.items():
        f.write(f"{k}={v}\n")
    f.write("\n# All datasets provided by current jarvis-tools:\n")
    for name in sorted(info.keys()):
        f.write(f"  - {name}\n")

(reports / "database_source.json").write_text(
    json.dumps(record, indent=2, ensure_ascii=False), encoding="utf-8"
)

print("\nWrote:", out)
print("Available datasets count:", len(info))
print("All datasets:", sorted(info.keys()))
