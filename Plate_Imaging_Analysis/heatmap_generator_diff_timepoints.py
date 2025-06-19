#!/usr/bin/env python
"""
Halo-clearing intensity heat-map

• reads  : series2_halo_clearing_cleaned_data.xls
• shows  : strains (rows) x 4 average-z conditions x 2(for both 48h and 72h)(columns)
• writes : halo_heatmap.png
"""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap  

# 1. load the cleaned data
df = pd.read_excel("series2_halo_clearing_cleaned_data.xlsx")


# 2. choose the four average-z columns
#    (untransformed on the left, transformed on the right)
cols = ["8h_pRLK120_25.0_avg_norm", "24h_pRLK120_25.0_avg_norm", "48h_pRLK120_25.0_avg_norm", "72h_pRLK120_25.0_avg_norm", "96h_pRLK120_25.0_avg_norm", 
        "8h_pRLK120_12.5_avg_norm","24h_pRLK120_12.5_avg_norm", "48h_pRLK120_12.5_avg_norm", "72h_pRLK120_12.5_avg_norm", "96h_pRLK120_12.5_avg_norm",
        "8h_pRLK139_25.0_avg_norm","24h_pRLK139_25.0_avg_norm", "48h_pRLK139_25.0_avg_norm",  "72h_pRLK139_25.0_avg_norm", "96h_pRLK139_25.0_avg_norm",
        "8h_pRLK139_12.5_avg_norm","24h_pRLK139_12.5_avg_norm", "48h_pRLK139_12.5_avg_norm", "72h_pRLK139_12.5_avg_norm", "96h_pRLK139_12.5_avg_norm"
        ]
cols = [c for c in cols if c in df.columns]  

heat = df.set_index("strain")[cols]

cmap = LinearSegmentedColormap.from_list(
    "sharp_bleach_red",
    ["#0075ff",  "#ffffff",  "#ff0000"],  
    N=25.6
)
norm = TwoSlopeNorm(
    vmin=heat.values.min(),
    vcenter=0,
    vmax=heat.values.max()
)

# 3. build heat-map
fig, ax = plt.subplots(figsize=(6, 12))
im = ax.imshow(heat.values, aspect="auto", cmap=cmap, norm=norm) 

# row / column labels
ax.set_yticks(range(len(heat.index)))
ax.set_yticklabels(heat.index)
ax.set_xticks(range(len(cols)))
ax.set_xticklabels(cols, rotation=45, ha="right")

# colour bar
cbar = plt.colorbar(im, ax=ax)
cbar.set_label("Average Intensity/Baseline Intensity")

ax.set_title("Halo-clearing average intensity per Strain at different timepoints")
plt.tight_layout()

# 4. save & display
plt.savefig("series2_diff_timepoints_heatmap.png", dpi=300)
plt.show()
print("✔ Saved series2_diff_timepoints_heatmap.png")
