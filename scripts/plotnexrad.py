import os
import re
import requests
import datetime
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pyart  # ensure pyart is installed in your environment

# ------------------------------
# Configuration
# ------------------------------
S3_BASE_URL = "https://noaa-nexrad-level2.s3.amazonaws.com/2025/08/28/KMOB/"
OUTPUT_DIR = "dirtonexradimage"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Map extent (same as tempmap.py)
# Format: [lon_min, lon_max, lat_min, lat_max]
MAP_EXTENT = [-87.5, -86.5, 30.5, 31.5]  # replace with your desired bounding box

# ------------------------------
# Fetch latest KMOB file
# ------------------------------
def get_latest_kmob_file(s3_base_url=S3_BASE_URL):
    response = requests.get(s3_base_url)
    response.raise_for_status()
    files = re.findall(r'KMOB\d+_\d+_V06', response.text)
    if not files:
        raise RuntimeError("No KMOB files found in S3 bucket")
    latest_file = sorted(files)[-1]
    file_url = f"{s3_base_url}{latest_file}"
    return file_url

try:
    nexrad_file_url = get_latest_kmob_file()
    print(f"Latest KMOB file URL: {nexrad_file_url}")
except Exception as e:
    print(f"Error fetching latest KMOB file: {e}")
    raise

# ------------------------------
# Download the file locally
# ------------------------------
local_filename = os.path.join(OUTPUT_DIR, nexrad_file_url.split("/")[-1])
if not os.path.exists(local_filename):
    print(f"Downloading {nexrad_file_url} ...")
    r = requests.get(nexrad_file_url, stream=True)
    with open(local_filename, 'wb') as f:
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)
else:
    print(f"File already exists locally: {local_filename}")

# ------------------------------
# Read and plot radar data
# ------------------------------
radar = pyart.io.read(local_filename)
fig = plt.figure(figsize=(10, 8))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent(MAP_EXTENT, crs=ccrs.PlateCarree())

# Add base map features
ax.add_feature(cfeature.STATES)
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS)

# Plot reflectivity
display = pyart.graph.RadarMapDisplay(radar)
display.plot_ppi_map(
    'reflectivity',
    0,
    ax=ax,
    title=f"KMOB Reflectivity {datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M UTC')}",
    cmap='pyart_NWSRef',
    colorbar_label='dBZ',
    vmin=-20,
    vmax=80
)

# Save the figure
output_file = os.path.join(OUTPUT_DIR, "nexradimage.png")
plt.savefig(output_file, dpi=150, bbox_inches='tight')
plt.close(fig)
print(f"Radar image saved to {output_file}")
