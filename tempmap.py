#!/usr/bin/env python3
"""
Florida Panhandle Temperature Map (NWS stations)
- Handles overlapping stations with slight offset
- Clamps labels to map boundaries
- Dynamic font scaling based on number of stations
- Output: Panhandle-Temp.png (1280x720)
"""

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
OUTFILE = "Panhandle-Temp.png"
FIGSIZE = (19.2, 10.8)  # inches
DPI = 100

print("Grabbing Stations")

HEADERS = {
    "User-Agent": "Milton Weather Dashboard (contact: your-email@domain.com)",
    "Accept": "application/geo+json",
}

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


print("Stations Successfully Grabbed")
print("Calculating Offsets")
# ----------------------------
# Auto-offset overlapping stations
# ----------------------------
def offset_points(points, min_distance=0.03):
    """
    Shift stations slightly if they are too close to avoid label overlap
    min_distance in degrees (~3 km)
    """
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
print("Calculating Font Sizes")
def choose_font_sizes(n_points: int):
    """Adjust font size based on number of stations"""
    if n_points >= 12:
        temp_fs = 28
        name_base = 10
    elif n_points >= 8:
        temp_fs = 34
        name_base = 12
    else:
        temp_fs = 40
        name_base = 14
    return temp_fs, name_base

def name_font_for(text: str, base: int):
    """Slightly shrink font for long names"""
    over = max(0, len(text) - 16)
    fs = max(10, int(base - 0.25 * over))
    return fs

# ----------------------------
# Plotting
# ----------------------------
def clamp(val, minval, maxval):
    return max(minval, min(maxval, val))


def plot_map(points):
    if not points:
        print("[INFO] No points to plot.")
        return

    points = offset_points(points)
    temp_fs, name_base = choose_font_sizes(len(points))

    fig = plt.figure(figsize=FIGSIZE, dpi=DPI)
    ax = plt.axes(projection=ccrs.Mercator())
    ax.set_extent(EXTENT, crs=ccrs.PlateCarree())

    # Base map
    ax.add_feature(cfeature.LAND.with_scale("10m"), facecolor="#d9d9d9")
    ax.add_feature(cfeature.OCEAN.with_scale("10m"), facecolor="#bcd7ff")
    ax.add_feature(cfeature.COASTLINE.with_scale("10m"), linewidth=0.8)
    ax.add_feature(cfeature.BORDERS.with_scale("10m"), linewidth=0.6)
    ax.add_feature(cfeature.STATES.with_scale("10m"), linewidth=0.6)

        # ----------------------------
    # Watermark
    # ----------------------------
    import matplotlib.image as mpimg
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    try:
        watermark = mpimg.imread("MiltonWxLogo.png")  # Ensure file exists in script folder
        # Create inset axes in bottom-right corner
        axins = inset_axes(ax,
                           width="25%",   # % of parent axes
                           height="25%",
                           loc='lower right',
                           borderpad=0.5)
        axins.imshow(watermark)
        axins.axis('off')  # hide axes ticks/frame
    except FileNotFoundError:
        print("[WARN] MiltonWxLogo.png not found. Skipping watermark.")

    outline = [pe.withStroke(linewidth=3.5, foreground="black", alpha=0.9)]

    for p in points:
        # Clamp points inside map
        x = clamp(p["lon"], EXTENT[0]+0.01, EXTENT[1]-0.01)
        y = clamp(p["lat"], EXTENT[2]+0.01, EXTENT[3]-0.01)

        # Optional dot marker
        ax.plot(x, y, marker="o", markersize=5,
                markerfacecolor="white", markeredgecolor="black",
                transform=ccrs.PlateCarree(), zorder=6)

        # Temperature label
        ax.text(
            x, y, f"{p['temp_f']}°",
            transform=ccrs.PlateCarree(),
            ha="center", va="bottom",
            fontsize=temp_fs, weight="bold", color="white",
            path_effects=[pe.withStroke(linewidth=3.5, foreground="black", alpha=0.9)],
            zorder=7,
        )

        # Station name
        name_fs = name_font_for(p["name"], name_base)
        ax.text(
            x, y - 0.025,
            p["name"],
            transform=ccrs.PlateCarree(),
            ha="center", va="top",
            fontsize=name_fs, color="white",
            path_effects=outline,
            zorder=7,
        )


    plt.savefig(OUTFILE, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved map as {OUTFILE}")

# ----------------------------
# Main
# ----------------------------
if __name__ == "__main__":
    points = load_all_stations(STATIONS)
    plot_map(points)
