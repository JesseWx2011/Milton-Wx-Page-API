"""
lightning.py

Fetches Blitzortung lightning data, plots on Cartopy using a local 10m countries shapefile,
and saves an image every 15 minutes.

Requires:
    pip install requests matplotlib cartopy pillow numpy
"""

import io
import re
import os
import tempfile
import time
import requests
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from PIL import Image
import numpy as np
from datetime import datetime, timezone
from zoneinfo import ZoneInfo

# ---- Configuration ----
PLACEFILE_URL = "https://saratoga-weather.org/USA-blitzortung/placefile.txt"
OUTPUT_FILE = "blitzortung_map.png"

# Path to local countries shapefile in repo root
REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
COUNTRIES_SHP = os.path.join(REPO_ROOT, "ne_10m_admin_0_countries.shp")

USE_ICON_IMAGE = False
ICON_INDEX = 9
SPRITE_TILE_SIZE = (30, 30)
MARKER_COLOR = "#ffff00"
MARKER_SIZE = 6
MARKER_ALPHA = 0.9
UPDATE_INTERVAL = 15 * 60  # 15 minutes in seconds

# ---- Helpers ----
def fetch_text(url):
    temp_path = os.path.join(tempfile.gettempdir(), "blitzortung_placefile.txt")
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
                      "AppleWebKit/537.36 (KHTML, like Gecko) "
                      "Chrome/140.0.0.0 Safari/537.36"
    }
    try:
        r = requests.get(url, headers=headers, timeout=10)
        r.raise_for_status()
        text = r.text
        with open(temp_path, "w", encoding="utf-8") as f:
            f.write(text)
        print(f"Downloaded placefile → {temp_path}")
        return text
    except Exception as e:
        print("Live download failed:", e)
        if os.path.exists(temp_path):
            print("Using cached placefile from temp folder.")
            with open(temp_path, "r", encoding="utf-8") as f:
                return f.read()
        else:
            raise RuntimeError("Could not fetch placefile and no cached copy available.")

def parse_placefile(text):
    lines = text.splitlines()
    meta = {"title": None, "refresh": None, "iconfile": None, "icons": []}
    icon_re = re.compile(r"^\s*Icon:\s*(.+)$", re.IGNORECASE)
    iconfile_re = re.compile(r"^\s*IconFile:\s*(.+)$", re.IGNORECASE)
    title_re = re.compile(r"^\s*Title:\s*(.+)$", re.IGNORECASE)
    refresh_re = re.compile(r"^\s*RefreshSeconds:\s*(\d+)", re.IGNORECASE)

    for ln in lines:
        ln = ln.strip()
        if not ln or ln.startswith(";"):
            continue
        m = iconfile_re.match(ln)
        if m:
            meta["iconfile"] = m.group(1).strip()
            continue
        m = title_re.match(ln)
        if m:
            meta["title"] = m.group(1).strip()
            continue
        m = refresh_re.match(ln)
        if m:
            meta["refresh"] = int(m.group(1))
            continue
        m = icon_re.match(ln)
        if m:
            rest = m.group(1).strip()
            parts = [p.strip() for p in rest.split(",", 5)]
            if len(parts) < 2:
                continue
            try:
                lat = float(parts[0])
                lon = float(parts[1])
            except ValueError:
                try:
                    lat = float(parts[1])
                    lon = float(parts[0])
                except Exception:
                    continue
            label = parts[-1] if len(parts) >= 3 else ""
            extra = parts[2:-1] if len(parts) > 3 else parts[2:-1]
            meta["icons"].append({"lat": lat, "lon": lon, "extra": extra, "label": label})
    return meta

def get_sprite_image(iconfile_url):
    try:
        resp = requests.get(iconfile_url, timeout=10)
        resp.raise_for_status()
        img = Image.open(io.BytesIO(resp.content)).convert("RGBA")
        return img
    except Exception as e:
        print("Could not fetch icon sprite:", e)
        return None

def extract_sprite_tile(sprite_image, index, tile_size):
    w, h = tile_size
    sheet_w, sheet_h = sprite_image.size
    tiles_per_row = sheet_w // w
    row = index // tiles_per_row
    col = index % tiles_per_row
    return sprite_image.crop((col*w, row*h, col*w+w, row*h+h))

# ---- Main plotting function ----
def plot_and_save():
    text = fetch_text(PLACEFILE_URL)
    meta = parse_placefile(text)
    icons = meta["icons"]
    if not icons:
        print("No lightning icons found.")
        return

    lats = [p["lat"] for p in icons]
    lons = [p["lon"] for p in icons]

    proj = ccrs.PlateCarree()
    fig = plt.figure(figsize=(14,10))
    ax = fig.add_axes([0.01, 0.05, 0.98, 0.92], projection=proj)

    # ---- Local 10m countries shapefile ----
    if os.path.exists(COUNTRIES_SHP):
        countries_reader = shpreader.Reader(COUNTRIES_SHP)
        ax.add_feature(
            cfeature.ShapelyFeature(
                countries_reader.geometries(),
                ccrs.PlateCarree(),
                facecolor="lightgreen",
                edgecolor="black"
            )
        )
    else:
        print(f"Warning: {COUNTRIES_SHP} not found. Using default Cartopy LAND + BORDERS.")
        ax.add_feature(cfeature.LAND)
        ax.add_feature(cfeature.BORDERS)

    # Optional ocean background
    ax.add_feature(cfeature.OCEAN, facecolor="lightblue")

    # Example bounds (adjust as needed)
    extent = [-126.329, -61.106299, 16.983, 61.106299]
    ax.set_extent(extent, crs=proj)

    # Title with CDT / UTC
    now_utc = datetime.now(timezone.utc)
    now_cdt = now_utc.astimezone(ZoneInfo("America/Chicago"))
    utc_str = now_utc.strftime("%Y-%m-%d %H:%M UTC")
    cdt_str = now_cdt.strftime("%Y-%m-%d %I:%M %p CDT")
    title = (meta.get("title") or "Blitzortung") + f"\n{cdt_str} / {utc_str}"
    ax.set_title(title, fontsize=14)

    # Plot points
    if USE_ICON_IMAGE and meta.get("iconfile"):
        sprite = get_sprite_image(meta["iconfile"])
        if sprite:
            tile = extract_sprite_tile(sprite, ICON_INDEX, SPRITE_TILE_SIZE)
            scale = 0.5
            tile = tile.resize((int(tile.size[0]*scale), int(tile.size[1]*scale)), Image.ANTIALIAS)
            arr_img = np.asarray(tile)
            for pt in icons:
                x, y = pt["lon"], pt["lat"]
                im = OffsetImage(arr_img, zoom=1)
                ab = AnnotationBbox(im, (x, y), frameon=False,
                                    xycoords=proj._as_mpl_transform(ax),
                                    boxcoords="offset points",
                                    pad=0)
                ax.add_artist(ab)
        else:
            ax.scatter(lons, lats, s=MARKER_SIZE, marker='o',
                       transform=proj, color=MARKER_COLOR, alpha=MARKER_ALPHA, zorder=3)
    else:
        ax.scatter(lons, lats, s=MARKER_SIZE, marker='o',
                   transform=proj, color=MARKER_COLOR, alpha=MARKER_ALPHA, zorder=3)

    # Save figure
    plt.savefig(OUTPUT_FILE, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved lightning map to {OUTPUT_FILE}")

# ---- Main loop ----
if __name__ == "__main__":
    while True:
        try:
            plot_and_save()
        except Exception as e:
            print("Error:", e)
        print(f"Waiting {UPDATE_INTERVAL/60:.0f} minutes for next update...")
        time.sleep(UPDATE_INTERVAL)
