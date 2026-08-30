import requests, json
H = {"User-Agent": "Mozilla/5.0 (research)"}
BASE = "https://jarvis.nist.gov/optimade/jarvisdft/v1/structures/"
r = requests.get(BASE, params={"filter": 'id STARTS WITH "dft_3d"', "page_limit": 2}, timeout=60, headers=H)
print("status", r.status_code)
d = r.json()
data = d.get("data", [])
print("returned", len(data))
for e in data[:2]:
    print("id:", e.get("id"))
    attrs = e.get("attributes", {})
    keys = sorted(attrs.keys())
    print("n_keys", len(keys))
    print("keys:", keys)
    # print kl-related
    for k in keys:
        if any(t in k.lower() for t in ["kl", "kappa", "therm", "conduct", "seebeck", "power", "cond", "mass", "gap", "phonon", "debye", "gruneisen", "eps", "slme"]):
            print("   PROP", k, "=", attrs[k])
