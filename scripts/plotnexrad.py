import os
import datetime
import requests
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs

from metpy.io import Level3File
from metpy.plots import add_metpy_logo, add_timestamp, colortables, USCOUNTIES
from metpy.calc import azimuth_range_to_lat_lon
from metpy.units import units

# Create output directory
output_dir = "dirtonexradimage"
os.makedirs(output_dir, exist_ok=True)

# S3 bucket URL for KMOB Level-III file (replace with latest URL logic if needed)
KMOB_file = "https://unidata-nexrad-level3.s3.amazonaws.com/MOB_N0B_2025_08_28_15_52_10"

# Download Level-III file
local_file = os.path.join(output_dir, os.path.basename(KMOB_file))
if not os.path.exists(local_file):
    r = requests.get(KMOB_file)
    r.raise_for_status()
    with open(local_file, "wb") as f:
        f.write(r.content)

# Open Level-III file with MetPy
f = Level3File(local_file)

# Grab first reflectivity product as an example
datadict = f.sym_block[0][0]
data = f.map_data(datadict['data'])

# Get azimuths and ranges
az = units.Quantity(np.array(datadict['start_az'] + [datadict['end_az'][-1]]), 'degrees')
rng = units.Quantity(np.linspace(0, f.max_range, data.shape[-1] + 1), 'kilometers')

# Central latitude and longitude
cent_lon, cent_lat = f.lon, f.lat

# Convert azimuth/range to lat/lon
xlocs, ylocs = azimuth_range_to_lat_lon(az, rng, cent_lon, cent_lat)

# Set up plot
fig = plt.figure(figsize=(10, 10))
spec = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(spec[0], projection=ccrs.LambertConformal())
ax.add_feature(USCOUNTIES, linewidth=0.5)

# Use MetPy colortable for reflectivity
norm, cmap = colortables.get_with_steps('NWSStormClearReflectivity', -20, 0.5)
ax.pcolormesh(xlocs, ylocs, data, norm=norm, cmap=cmap, transform=ccrs.PlateCarree())

# Focus on your area (example: Milton, FL)
ax.set_extent([-88.38, -86.337, 30.7, 30.182])
ax.set_aspect('equal', 'datalim')

# Add MetPy logo and timestamp
add_metpy_logo(fig, 190, 85, size='large')
add_timestamp(ax, f.metadata['prod_time'], y=0.02, high_contrast=True)

# Save image
out_file = os.path.join(output_dir, "nexradimage.png")
plt.savefig(out_file, dpi=150, bbox_inches='tight')
plt.close()

print(f"Radar image saved to: {out_file}")
