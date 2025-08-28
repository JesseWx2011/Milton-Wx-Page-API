#!/usr/bin/env python3
"""
Fetch latest KMOB NEXRAD file (Level-III or Level-II .gz/.ar2v) and plot reflectivity
cropped to Florida Panhandle extent (same as tempmap.py).
Saves output to docs/radar_images/KMOB_radar_<timestamp>.png
"""

import os
import re
import sys
import time
import tempfile
import datetime
from datetime import timezone
import requests
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib import patheffects as pe
import matplotlib.image as mpimg
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Optional pyart support for Level-II
try:
    import pyart
    HAVE_PYART = True
except Exception:
    HAVE_PYART = False
    print("[warn] pyart not available: Level-II (.ar2v/.gz) support disabled. Install pyart if you need Level-II support.")

# ----------------------------
# User config
# ----------------------------
RADAR_STN = "KMOB"
# Panhandle extent (lon_min, lon_max, lat_min, lat_max)
EXTENT = [-88.4, -86.1, 30.3, 31.1]
OUTDIR = "docs/radar_images"
os.makedirs(OUTDIR, exist_ok=True)

# How many days to search (0 = today)
SEARCH_DAYS_BACK = 2

# If automatic discovery fails, set this to a known-good URL (optional)
FALLBACK_URL = None
# Example: FALLBACK_URL = "https://noaa-nexrad-level3.s3.amazonaws.com/2025/08/27/KMOB/KMOB20250827_2110.gz"

# Tile / plotting options
VMIN, VMAX = -32, 84
SHADING = 'auto'  # 'auto', 'flat', 'nearest'
DPI = 150

# ----------------------------
# Helpers: discover latest Level-III (multiple strategies)
# ----------------------------
def list_ncei_level3_files(year, month, day, station):
    """
    Try NCEI 'access' directory listing for Level-3:
    https://www.ncei.noaa.gov/data/nexrad-level-3/access/YYYY/MM/DD/STATION/
    Returns list of full file URLs (may be [] if not accessible).
    """
    base = f"https://www.ncei.noaa.gov/data/nexrad-level-3/access/{year}/{month:02d}/{day:02d}/{station}/"
    try:
        r = requests.get(base, timeout=15)
        r.raise_for_status()
        html = r.text
        # look for href= links to .gz or .nids or similar
        links = re.findall(r'href=[\'"]([^\'"]+\.(?:gz|ar2v|nids|nids.ascii|nc))', html, flags=re.IGNORECASE)
        urls = []
        for link in links:
            if link.startswith("http"):
                urls.append(link)
            else:
                urls.append(base + link)
        return urls
    except Exception as e:
        # not available / listing blocked
        return []

def list_s3_level3_files(year, month, day, station):
    """
    Try the common S3 bucket listing that some NOAA mirrors use.
    Pattern: https://noaa-nexrad-level3.s3.amazonaws.com/YYYY/MM/DD/STATION/
    Not guaranteed to be enabled, but attempt it.
    """
    base = f"https://noaa-nexrad-level3.s3.amazonaws.com/{year}/{month:02d}/{day:02d}/{station}/"
    try:
        r = requests.get(base, timeout=15)
        r.raise_for_status()
        html = r.text
        links = re.findall(r'href=[\'"]([^\'"]+\.(?:gz|ar2v|nids|nc))', html, flags=re.IGNORECASE)
        urls = []
        for link in links:
            if link.startswith("http"):
                urls.append(link)
            else:
                urls.append(base + link)
        return urls
    except Exception:
        return []

def discover_latest_kmob():
    """Try several sources to find the latest KMOB Level-III or Level-II file."""
    now_utc = datetime.datetime.utcnow()
    candidates = []
    for dback in range(0, SEARCH_DAYS_BACK + 1):
        day = now_utc - datetime.timedelta(days=dback)
        y, m, daynum = day.year, day.month, day.day
        # 1) NCEI access listing
        urls = list_ncei_level3_files(y, m, daynum, RADAR_STN)
        if urls:
            candidates.extend(urls)
        # 2) noaa-nexrad-level3 S3
        urls2 = list_s3_level3_files(y, m, daynum, RADAR_STN)
        if urls2:
            candidates.extend(urls2)

    # Keep unique and sort (newest by filename/time heuristic)
    candidates = list(dict.fromkeys(candidates))
    if not candidates:
        return None

    # Heuristic: filenames often contain timestamps. We'll pick the lexicographically largest .gz/.nc
    candidates_sorted = sorted(candidates, reverse=True)
    # Prefer .gz /.ar2v /.nc in that order
    for ext in ('.gz', '.ar2v', '.nc', '.nids', '.nids.ascii'):
        for c in candidates_sorted:
            if c.lower().endswith(ext):
                return c
    return candidates_sorted[0]

# ----------------------------
# Read reflectivity from Level-III xarray dataset
# ----------------------------
def get_reflectivity_level3_from_ds(ds):
    # Try a few possible variable/coordinate names
    lat = lon = None
    for name in ('latitude','radar_latitude','radar_lat','site_latitude'):
        if name in ds:
            lat = float(ds[name].values.ravel()[0])
            break
    for name in ('longitude','radar_longitude','radar_lon','site_longitude'):
        if name in ds:
            lon = float(ds[name].values.ravel()[0])
            break

    # Try reflectivity vars
    if 'BaseReflectivityDR' in ds:
        refl = ds['BaseReflectivityDR'].values
    elif 'BaseReflectivityDR_RAW' in ds:
        raw = ds['BaseReflectivityDR_RAW'].values
        a = np.asarray(raw, dtype=np.uint8)
        nodata_mask = (a == 0)
        val = a & 0x7F
        refl = val.astype(np.float32) - 32.0
        refl[nodata_mask] = np.nan
    elif 'reflectivity' in ds:
        refl = ds['reflectivity'].values
    else:
        raise KeyError("No recognized reflectivity variable found in Level-III dataset")

    # Azimuths / gates
    az = np.asarray(ds['azimuth'].values) if 'azimuth' in ds else np.linspace(0, 360, refl.shape[0], endpoint=False)
    gate = np.asarray(ds['gate'].values) if 'gate' in ds else np.arange(refl.shape[1]) * 250.0

    # Try refTime timestamp
    scan_time = None
    if 'refTime' in ds:
        try:
            scan_time = xr.conventions.decode_cf_datetime(ds['refTime'].values, units='s')
        except Exception:
            try:
                scan_time = np.datetime64(ds['refTime'].values).astype('datetime64[ms]').astype('datetime64[ms]').tolist()
            except Exception:
                scan_time = None

    return refl, az, gate, lat, lon, scan_time

# ----------------------------
# Read Level-II using Py-ART (if available)
# ----------------------------
def get_reflectivity_level2_with_pyart(file_path):
    if not HAVE_PYART:
        raise RuntimeError("pyart not installed; cannot parse Level-II file")

    # allow http(s) gz download
    tmp_path = None
    try:
        if file_path.startswith("http") and file_path.endswith(".gz"):
            r = requests.get(file_path, stream=True, timeout=30)
            r.raise_for_status()
            fd, tmp_path = tempfile.mkstemp(suffix=".gz")
            with os.fdopen(fd, "wb") as f:
                f.write(r.content)
            file_to_read = tmp_path
        else:
            file_to_read = file_path

        radar = pyart.io.read(file_to_read)
        # pick reflectivity field
        if 'reflectivity' in radar.fields:
            refl = radar.fields['reflectivity']['data']
        elif 'reflectivity_horizontal' in radar.fields:
            refl = radar.fields['reflectivity_horizontal']['data']
        else:
            raise KeyError("No reflectivity field found in Level-II read")

        az = radar.azimuth['data']
        gate = radar.range['data']
        lat0 = float(radar.latitude['data'][0])
        lon0 = float(radar.longitude['data'][0])
        # time parsing:
        try:
            scan_time = None
            if hasattr(radar, 'time'):
                scan_time = None  # pyart time parsing can be added here if needed
        except Exception:
            scan_time = None
        return refl, az, gate, lat0, lon0, scan_time
    finally:
        if tmp_path:
            try:
                os.remove(tmp_path)
            except Exception:
                pass

# ----------------------------
# Compute lat/lon grid (simple geodesic option)
# ----------------------------
def compute_latlon_grid(lat0, lon0, azimuths_deg, gates_m):
    # vectorized forward using simple spherical approx (good for small ranges)
    # Use pyproj.Geod if available for better accuracy
    try:
        from pyproj import Geod
        geod = Geod(ellps='WGS84')
        az_rad = np.deg2rad(azimuths_deg)
        # create 2D arrays
        az2d = np.repeat(azimuths_deg[:, np.newaxis], len(gates_m), axis=1)
        gate2d = np.repeat(gates_m[np.newaxis, :], len(azimuths_deg), axis=0)
        lon2d = np.zeros_like(az2d, dtype=float)
        lat2d = np.zeros_like(az2d, dtype=float)
        # vectorized loop (pyproj does not vectorize directly)
        for i, az in enumerate(azimuths_deg):
            lon_arr, lat_arr, _ = geod.fwd(np.full_like(gates_m, lon0),
                                           np.full_like(gates_m, lat0),
                                           np.full_like(gates_m, az),
                                           gate2d[i])
            lon2d[i, :] = lon_arr
            lat2d[i, :] = lat_arr
        return lon2d, lat2d
    except Exception:
        # fallback: approximate using simple haversine forward (less accurate)
        az2d = np.repeat(azimuths_deg[:, np.newaxis], len(gates_m), axis=1)
        gate2d = np.repeat(gates_m[np.newaxis, :], len(azimuths_deg), axis=0)
        R = 6371000.0
        lat0_rad = np.deg2rad(lat0)
        lon0_rad = np.deg2rad(lon0)
        lat2d = np.arcsin(np.sin(lat0_rad) * np.cos(gate2d / R) +
                          np.cos(lat0_rad) * np.sin(gate2d / R) * np.cos(np.deg2rad(az2d)))
        lon2d = lon0_rad + np.arctan2(np.sin(np.deg2rad(az2d)) * np.sin(gate2d / R) * np.cos(lat0_rad),
                                      np.cos(gate2d / R) - np.sin(lat0_rad) * np.sin(lat2d))
        return np.rad2deg(lon2d), np.rad2deg(lat2d)

# ----------------------------
# Plot function: crop to EXTENT
# ----------------------------
def plot_reflectivity(lon2d, lat2d, refl2d, radar_lon, radar_lat, scan_time=None):
    fig = plt.figure(figsize=(12.8, 7.2), dpi=DPI)  # 1280x720
    ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree())
    ax.set_extent(EXTENT, crs=ccrs.PlateCarree())

    # base layers
    ax.add_feature(cfeature.LAND.with_scale("10m"), facecolor="#d9d9d9")
    ax.add_feature(cfeature.OCEAN.with_scale("10m"), facecolor="#bcd7ff")
    ax.add_feature(cfeature.COASTLINE.with_scale("10m"), linewidth=0.6)
    ax.add_feature(cfeature.STATES.with_scale("10m"), linewidth=0.5)

    # colormap: use matplotlib's default radar-like palette
    cmap = plt.get_cmap("turbo")
    norm = plt.Normalize(vmin=VMIN, vmax=VMAX)

    # mask invalid
    data = np.where(np.isfinite(refl2d), refl2d, np.nan)
    mesh = ax.pcolormesh(lon2d, lat2d, data, transform=ccrs.PlateCarree(),
                         cmap=cmap, norm=norm, shading=SHADING, zorder=4)

    # radar marker
    ax.plot(radar_lon, radar_lat, marker='^', color='k', markersize=8, transform=ccrs.PlateCarree(), zorder=6)

    cb = plt.colorbar(mesh, ax=ax, orientation='vertical', pad=0.02, fraction=0.046)
    cb.set_label("Reflectivity (dBZ)")
    cb.ax.tick_params(labelsize=8)

    # watermark
    try:
        watermark = mpimg.imread("MiltonWxLogo.png")
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        axins = inset_axes(ax, width="18%", height="18%", loc='lower right', borderpad=0.5)
        axins.imshow(watermark)
        axins.axis('off')
    except FileNotFoundError:
        pass

    # timestamp
    if scan_time is not None:
        if isinstance(scan_time, (np.datetime64, datetime.datetime)):
            try:
                dt_utc = np.datetime64(scan_time).astype('datetime64[s]').astype(datetime.datetime)
            except Exception:
                dt_utc = scan_time
            utc = dt_utc if isinstance(dt_utc, datetime.datetime) else None
            if utc is not None:
                utc = utc.replace(tzinfo=timezone.utc)
                local = utc.astimezone(datetime.timezone(datetime.timedelta(hours=-6)))  # rough US Central fallback
                label = f"Scan: {local.strftime('%Y-%m-%d %I:%M:%S %p %Z')} / {utc.strftime('%Y-%m-%d %H:%M UTC')}"
            else:
                label = f"Scan: {scan_time}"
        else:
            label = f"Scan: {scan_time}"
        ax.text(0.01, 0.01, label, transform=ax.transAxes, fontsize=9,
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

    # save
    ts = datetime.datetime.utcnow().strftime("%Y%m%d_%H%M%S")
    outname = os.path.join(OUTDIR, f"{RADAR_STN}_reflectivity_{ts}.png")
    plt.tight_layout()
    plt.savefig(outname, dpi=DPI)
    plt.close(fig)
    print(f"[info] saved {outname}")

# ----------------------------
# Main flow
# ----------------------------
def main():
    print("[info] discovering latest KMOB file...")
    url = discover_latest_kmob()
    if not url:
        if FALLBACK_URL:
            url = FALLBACK_URL
            print(f"[warn] discovery failed, using fallback URL: {url}")
        else:
            print("[error] could not find a KMOB file automatically. Set FALLBACK_URL to a valid file URL and retry.")
            sys.exit(1)

    print(f"[info] using {url}")

    # if remote .gz Level-II, try pyart (if installed), else try Level-III xarray open
    lower = url.lower()
    try:
        if (lower.endswith(".ar2v") or lower.endswith(".gz")) and HAVE_PYART:
            print("[info] attempting Level-II parse via pyart...")
            refl, az, gate, lat0, lon0, scan_time = get_reflectivity_level2_with_pyart(url)
        else:
            # Try opening with xarray; many Level-III datasets are netcdf/CF-accessible
            print("[info] attempting Level-III/CF parse via xarray...")
            ds = xr.open_dataset(url, decode_times=False)
            refl, az, gate, lat0, lon0, scan_time = get_reflectivity_level3_from_ds(ds)
    except Exception as e:
        print(f"[error] parsing radar file failed: {e}")
        # try fallback to pyart if possible
        if HAVE_PYART and (lower.endswith(".ar2v") or lower.endswith(".gz")):
            try:
                refl, az, gate, lat0, lon0, scan_time = get_reflectivity_level2_with_pyart(url)
            except Exception as ee:
                print(f"[error] pyart fallback failed: {ee}")
                sys.exit(1)
        else:
            sys.exit(1)

    # compute lat/lon grid
    try:
        gates_m = np.asarray(gate)
        azs = np.asarray(az)
    except Exception:
        # fallback sizes
        gates_m = np.arange(refl.shape[1]) * 250.0
        azs = np.linspace(0, 360, refl.shape[0], endpoint=False)

    lon2d, lat2d = compute_latlon_grid(lat0, lon0, azs, gates_m)

    # ensure shape matches (attempt to resize or slice)
    try:
        if lon2d.shape != refl.shape:
            # do a simple interpolation/resample if needed
            from scipy.ndimage import zoom as ndzoom
            zoom_factors = (lon2d.shape[0]/refl.shape[0], lon2d.shape[1]/refl.shape[1])
            refl2 = ndzoom(refl, zoom_factors, order=1)
        else:
            refl2 = refl
    except Exception:
        refl2 = refl

    # Finally plot cropped to EXTENT
    plot_reflectivity(lon2d, lat2d, refl2, lon0, lat0, scan_time=scan_time)

if __name__ == "__main__":
    main()
