import json
import yaml
import pandas as pd
import matplotlib.pyplot as plt
import os

# Example station data (replace with your API fetch later)
stations = [
    {"id": "KFLMILTO379", "name": "Milton, FL", "lat": 30.63, "lon": -87.04},
    {"id": "KMOB", "name": "Mobile, AL", "lat": 30.68, "lon": -88.24},
]

# Convert to DataFrame
df = pd.DataFrame(stations)

# Ensure output directories exist
os.makedirs("docs/data", exist_ok=True)
os.makedirs("docs/images", exist_ok=True)

# Write JSON
with open("docs/data/stations.json", "w") as f:
    json.dump(stations, f, indent=2)

# Write YAML
with open("docs/data/stations.yml", "w") as f:
    yaml.dump(stations, f, sort_keys=False)

# Plot stations (simple scatter plot for now)
plt.figure(figsize=(6, 6))
plt.scatter(df["lon"], df["lat"], marker="o", c="red", s=100)
for _, row in df.iterrows():
    plt.text(row["lon"] + 0.05, row["lat"] + 0.05, row["name"], fontsize=8)
plt.title("Weather Stations")
plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.grid(True)
plt.savefig("docs/images/stations.png", dpi=150)
plt.close()
