import matplotlib.pyplot as plt
import json
from datetime import datetime

# Load your weather data (example from the generated JSON)
with open("docs/data/stations.json", "r") as f:
    data = json.load(f)

# Just make a quick visualization
temps = [station["temperature"] for station in data["stations"]]
names = [station["id"] for station in data["stations"]]

plt.figure(figsize=(8, 5))
plt.bar(names, temps)
plt.title(f"Station Temperatures ({datetime.utcnow().strftime('%Y-%m-%d %H:%M UTC')})")
plt.ylabel("Temperature (°F)")

plt.tight_layout()
plt.savefig("docs/data/weather_chart.png")  # save into docs/ so GitHub Pages can serve it
