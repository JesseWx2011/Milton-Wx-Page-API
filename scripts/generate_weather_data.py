#!/usr/bin/env python3
"""
Fetch Ambient + NWS latest observations and write YAML/JSON into docs/data/
Intended to run in GitHub Actions (secrets provided as env vars).
"""

import os
import requests
import yaml
import json
from datetime import datetime, timezone

# --- Config ---
AMBIENT_APP_KEY = os.environ.get("AMBIENT_APP_KEY")
AMBIENT_API_KEY = os.environ.get("AMBIENT_API_KEY")

AMBIENT_URL = (
    f"https://api.ambientweather.net/v1/devices?"
    f"applicationKey={AMBIENT_APP_KEY}&apiKey={AMBIENT_API_KEY}"
)

# NWS stations to include
NWS_STATIONS = ["KPNS", "KMOB", "KCEW", "K54J", "KJKA", "KNPA"]

# Output files (must be under docs/ to be served by GitHub Pages)
OUT_YAML = "docs/data/stations.yml"
OUT_JSON = "docs/data/stations.json"

HEADERS = {
    "User-Agent": "Milton Weather Dashboard (GitHub Actions)",
    "Accept": "application/geo+json",
}
# -----------------------

def fetch_ambient():
    if not (AMBIENT_APP_KEY and AMBIENT_API_KEY):
        print("Ambient keys not set, skipping Ambient fetch.")
        return []
    try:
        r = requests.get(AMBIENT_URL, timeout=15)
        r.raise_for_status()
        devices = r.json()
    except Exception as e:
        print("Ambient fetch failed:", e)
        return []

    out = []
    for d in devices:
        try:
            last = d.get("lastData", {}) or {}
            info = d.get("info", {}) or {}
            coords = info.get("coords", {}).get("coords", {}) or {}
            lat = coords.get("lat")
            lon = coords.get("lon")
            tempf = last.get("tempf")
            name = info.get("name") or info.get("location") or d.get("macAddress")
            if lat is None or lon is None:
                continue
            out.append({
                "source": "ambient",
                "id": d.get("macAddress"),
                "name": name,
                "lat": float(lat),
                "lon": float(lon),
                "temp_f": float(tempf) if tempf is not None else None,
                "timestamp": last.get("date") or last.get("dateutc")
            })
        except Exception as e:
            print("Ambient parse error:", e)
    return out

def fetch_nws_station(station_id):
    url = f"https://api.weather.gov/stations/{station_id}/observations/latest"
    try:
        r = requests.get(url, headers=HEADERS, timeout=12)
        r.raise_for_status()
        data = r.json()
    except Exception as e:
        print(f"NWS fetch failed for {station_id}: {e}")
        return None

    props = data.get("properties", {}) or {}
    geom = data.get("geometry", {}) or {}
    temp_c = props.get("temperature", {}).get("value")
    if temp_c is None:
        return None
    try:
        temp_f = round(temp_c * 9/5 + 32, 1)
    except Exception:
        temp_f = None

    coords = geom.get("coordinates")
    if not coords or len(coords) < 2:
        return None
    lon, lat = float(coords[0]), float(coords[1])
    name = props.get("stationName") or props.get("stationId") or station_id
    ts = props.get("timestamp") or props.get("validTime") or None

    return {
        "source": "nws",
        "id": station_id,
        "name": name,
        "lat": lat,
        "lon": lon,
        "temp_f": temp_f,
        "timestamp": ts
    }

def main():
    all_points = []

    # Ambient
    ambient_points = fetch_ambient()
    all_points.extend(ambient_points)

    # NWS
    for sid in NWS_STATIONS:
        st = fetch_nws_station(sid)
        if st:
            all_points.append(st)

    # Add generation timestamp
    generated = datetime.now(timezone.utc).astimezone().isoformat()
    payload = {
        "generated_at": generated,
        "count": len(all_points),
        "stations": all_points
    }

    # Ensure docs/data exists
    os.makedirs(os.path.dirname(OUT_YAML), exist_ok=True)

    # Write YAML
    with open(OUT_YAML, "w", encoding="utf-8") as f:
        yaml.safe_dump(payload, f, sort_keys=False, allow_unicode=True)

    # Also write JSON (handy)
    with open(OUT_JSON, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)

    print(f"Wrote {OUT_YAML} and {OUT_JSON}")

if __name__ == "__main__":
    main()
