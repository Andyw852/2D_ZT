"""Phase D - Step 4: 验证快照完整性 + SHA256 + 基础校验。"""
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

root = Path(__file__).resolve().parents[1]
raw_dir = root / "data" / "raw" / "jarvis"
snapshot = raw_dir / "dft_2d_snapshot.json"

st = snapshot.stat()
size = st.st_size
mtime = datetime.fromtimestamp(st.st_mtime, tz=timezone.utc).isoformat()

h = hashlib.sha256()
with snapshot.open("rb") as f:
    for chunk in iter(lambda: f.read(1 << 20), b""):
        h.update(chunk)
sha = h.hexdigest()

# 写入 SHA256SUMS（与 sha256sum 输出格式一致）
sums_file = raw_dir / "SHA256SUMS.txt"
with sums_file.open("w") as f:
    f.write(f"{sha}  {snapshot.name}\n")

with snapshot.open(encoding="utf-8") as f:
    records = json.load(f)

jids = [r.get("_jarvis_jid") for r in records]
formulas = [r.get("_jarvis_formula") for r in records]
unique_jids = len(set(jids))
n_none_jid = sum(1 for j in jids if j is None)

report = {
    "execution_date": datetime.now(timezone.utc).isoformat(),
    "file": str(snapshot),
    "file_size_bytes": size,
    "file_size_mb": round(size / 1e6, 3),
    "mtime": mtime,
    "sha256": sha,
    "record_count": len(records),
    "unique_jids": unique_jids,
    "duplicate_jids": len(records) - unique_jids,
    "n_none_jid": n_none_jid,
    "sample_jids": jids[:5],
    "sample_formulas": formulas[:5],
    "type_is_list": isinstance(records, list),
}
out = raw_dir / "archive_validation.json"
out.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
print(json.dumps(report, indent=2, ensure_ascii=False))
