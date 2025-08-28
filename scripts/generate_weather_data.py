#!/usr/bin/env python3
"""
Generate Panhandle Temperature map and save to docs/images/Panhandle-Temp.png
Uses NWS stations (and you can add Ambient fetch if desired).
"""

import os
import requests
import matplotlib.pyplot as plt
from matplotlib import patheffects as pe
import matplotlib.image as mpimg
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime, timezone

# ----------------------------
# Config
# ----------------------------
STATIONS = ["KPNS", "KMOB", "KCEW", "K54J", "KJKA", "KNPA"]  # Add/remove as needed
EXTENT = [-88.4, -86.1, 30.3, 31.1]  # [lon_min, lon_max, lat_min, lat_max]
OUTFILE = "docs/images/Panhandle-Temp.png"   # <<--- save inside docs/images for GitHub Pages
FIGSIZE = (19.2, 10.8)  # inches -> 1920x1080 at DPI=100
DPI = 100

HEADERS = {
    "User-Agent": "Milton Weather Dashboard (GitHub Actions)",
    "Accept": "application/geo+json",
}

# Ensure output directories exist
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)
os.makedirs("docs/data", exist_ok=True)  # if you also write data later

# ----------------------------
# Fetch NWS Observations
# ----------------------------
def fetch_latest_obs(station_id: str):
    url = f"https://api.weather.gov/stations/{station_id}/observations/latest"
    try:
        r = requests.get(url, headers=HEADERS, timeout=12)
        r.raise_for_status()
        data = r.json()

        props = data.get("properties", {}) or {}
        geom = data.get("geometry", {}) or {}

        # Temperature in °C -> °F
        temp_c = props.get("temperature", {}).get("value")
        if temp_c is None:
            return None
        temp_f = round(temp_c * 9/5 + 32)

        # Coordinates at top-level geometry
        coords = geom.get("coordinates")
        if not coords or len(coords) < 2:
            return None
        lon, lat = float(coords[0]), float(coords[1])

        name = props.get("stationName") or props.get("stationId") or station_id

        return {"id": station_id, "name": name, "lon": lon, "lat": lat, "temp_f": temp_f}
    except Exception as e:
        print(f"[WARN] {station_id}: {e}")
        return None

def load_all_stations(station_ids):
    obs = []
    for sid in station_ids:
        item = fetch_latest_obs(sid)
        if item:
            obs.append(item)
    return obs

# ----------------------------
# Offset overlapping stations
# ----------------------------
def offset_points(points, min_distance=0.03):
    adjusted = []
    for i, p1 in enumerate(points):
        x, y = p1["lon"], p1["lat"]
        shift_x, shift_y = 0, 0
        for p2 in adjusted:
            dx = x + shift_x - p2["lon"]
            dy = y + shift_y - p2["lat"]
            dist = (dx**2 + dy**2)**0.5
            if dist < min_distance:
                shift_x += min_distance/2
                shift_y += min_distance/2
        adjusted.append({**p1, "lon": x + shift_x, "lat": y + shift_y})
    return adjusted

# ----------------------------
# Font utilities
# ----------------------------
def choose_font_sizes(n_points: int):
    if n_points >= 12:
        temp_fs = 36
        name_base = 12
    elif n_points >= 8:
        temp_fs = 42
        name_base = 14
    else:
        temp_fs = 50
        name_base = 16
    return temp_fs, name_base

def name_font_for(text: str, base: int):
    over = max(0, len(text) - 16)
    fs = max(10, int(base - 0.25 * over))
    return fs

def clamp(val, minval, maxval):
    return max(minval, min(maxval, val))

# ----------------------------
# Plot map
# ----------------------------
def plot_map(points):
    if not points:
        print("[INFO] No points to plot.")
        return

    points = offset_points(points)
    temp_fs, name_base = choose_font_sizes(len(points))

    fig = plt.figure(figsize=FIGSIZE, dpi=DPI)
    ax = plt.axes(projection=ccrs.Mercator())
    ax.set_extent(EXTENT, crs=ccrs.PlateCarree())

    # Base map (10m for crisp lines)
    ax.add_feature(cfeature.LAND.with_scale("10m"), facecolor="#d9d9d9")
    ax.add_feature(cfeature.OCEAN.with_scale("10m"), facecolor="#bcd7ff")
    ax.add_feature(cfeature.COASTLINE.with_scale("10m"), linewidth=0.8)
    ax.add_feature(cfeature.BORDERS.with_scale("10m"), linewidth=0.6)
    ax.add_feature(cfeature.STATES.with_scale("10m"), linewidth=0.6)

    # Watermark (bottom-right)
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    try:
        watermark = mpimg.imread("MiltonWxLogo.png")
        axins = inset_axes(ax, width="18%", height="18%", loc='lower right', borderpad=0.5)
        axins.imshow(watermark)
        axins.axis('off')
    except FileNotFoundError:
        print("[WARN] MiltonWxLogo.png not found. Skipping watermark.")

    outline = [pe.withStroke(linewidth=3.5, foreground="black", alpha=0.9)]

    for p in points:
        x = clamp(p["lon"], EXTENT[0]+0.01, EXTENT[1]-0.01)
        y = clamp(p["lat"], EXTENT[2]+0.01, EXTENT[3]-0.01)

        # Marker
        ax.plot(x, y, marker="o", markersize=6,
                markerfacecolor="white", markeredgecolor="black",
                transform=ccrs.PlateCarree(), zorder=6)

        # Temperature label
        ax.text(x, y, f"{p['temp_f']}°",
                transform=ccrs.PlateCarree(),
                ha="center", va="bottom",
                fontsize=temp_fs, weight="bold", color="white",
                path_effects=outline, zorder=7)

        # Station name
        name_fs = name_font_for(p["name"], name_base)
        ax.text(x, y - 0.025, p["name"],
                transform=ccrs.PlateCarree(),
                ha="center", va="top",
                fontsize=name_fs, color="white",
                path_effects=outline, zorder=7)

    # Centered title at top
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    plt.title(f"Current Temperatures — {now}", fontsize=18, weight="bold", loc="center", pad=20)

    # Save final image (in docs/images so GitHub Pages serves it)
    plt.savefig(OUTFILE, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved map as {OUTFILE}")

# ----------------------------
# Main
# ----------------------------
if __name__ == "__main__":
    points = load_all_stations(STATIONS)
    plot_map(points)
