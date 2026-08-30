"""下载 Materials Project elasticity (含热导率) + summary (带隙) 到 mp_kappaL/raw。"""
import json, time, urllib.request, urllib.parse
from pathlib import Path

KEY = "UxK6s0AbvNiBDzrzOqH5rgTtmwW4nc22"
BASE = "https://api.materialsproject.org"
root = Path(__file__).resolve().parents[1] / "mp_kappaL"
raw = root / "raw"
raw.mkdir(parents=True, exist_ok=True)

def get(path, params):
    url = BASE + path + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"X-API-KEY": KEY, "User-Agent": "research"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            print("  retry", attempt, type(e).__name__, str(e)[:120])
            time.sleep(2 * (attempt + 1))
    raise RuntimeError("download failed: " + url)

# ---- 1. elasticity ----
EL_FIELDS = ("material_id,formula_pretty,structure,nsites,elements,density,volume,"
             "bulk_modulus,shear_modulus,debye_temperature,thermal_conductivity,"
             "sound_velocity,universal_anisotropy,homogeneous_poisson")
all_el = []
skip = 0
first = get("/materials/elasticity/", {"_fields": "material_id", "_limit": 1})
total = first["meta"]["total_doc"]
print("elasticity total_doc:", total)
while len(all_el) < total:
    d = get("/materials/elasticity/", {"_fields": EL_FIELDS, "_limit": 1000, "_skip": skip})
    data = d.get("data", [])
    if not data:
        break
    all_el.extend(data)
    skip += 1000
    print("  elapsed elasticity:", len(all_el), "/", total)
    if len(data) < 1000:
        break

(raw / "elasticity_all.json").write_text(json.dumps(all_el), encoding="utf-8")
print("saved elasticity_all.json:", len(all_el))

# ---- 2. summary band gap ----
mids = [x["material_id"] for x in all_el]
print("material_ids:", len(mids))
summ = {}
batch = 500
for i in range(0, len(mids), batch):
    chunk = mids[i:i+batch]
    d = get("/materials/summary/", {
        "_fields": "material_id,formula_pretty,band_gap,is_metal,is_gap_direct,efermi",
        "material_ids": ",".join(chunk)})
    for x in d.get("data", []):
        summ[x["material_id"]] = x
    print("  elapsed summary:", len(summ), "/", len(mids))

(raw / "summary_bandgap.json").write_text(json.dumps(summ), encoding="utf-8")
print("saved summary_bandgap.json:", len(summ))
