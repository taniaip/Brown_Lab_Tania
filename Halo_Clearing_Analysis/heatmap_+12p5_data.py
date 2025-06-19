#!/usr/bin/env python3
# PETase halo assay heat-map  –  updated for 25 mM & 12.5 mM
# -----------------------------------------------------------
import re
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# 1. load data ------------------------------------------------------------
df_raw = (
    pd.read_csv("clreared_halo24+42_clearinng48.csv")
      .set_index("strain")
)

# 2. pick out colony-intensity columns by regex ---------------------------
#    Ha. → 25 mM      Cl. → 12.5 mM
pattern = re.compile(r"^(Ha|Cl)\.(X\d+)(b?)\.colony\.value\.(\d+)h$")
col_cols = [c for c in df_raw.columns if pattern.match(c)]
df_vals  = df_raw[col_cols].apply(pd.to_numeric, errors="coerce")

# 3. aggregate replicates:  conc × status × time -------------------------
#    build grouping dict
groups = {}
for col in df_vals.columns:
    src, _, b_flag, t = pattern.match(col).groups()
    conc   = "25mM"   if src == "Ha" else "12.5mM"
    status = "UnT" if b_flag else "T"
    key    = (conc, status, f"{t}h")
    groups.setdefault(key, []).append(col)

# define valid time‐points per concentration
times_by_conc = {
    "25mM": ["24h", "42h"],    # halo assay
    "12.5mM": ["48h"],         # clearing assay
}

concs  = ["25mM", "12.5mM"]
status = ["UnT", "T"]

# build agg dict of pandas.Series
agg = {}
for conc in concs:
    for st in status:
        for t in times_by_conc[conc]:
            cols = groups.get((conc, st, t), [])
            if cols:
                agg[f"{conc}_{st}_{t}"] = df_vals[cols].mean(axis=1, skipna=True)
            else:
                # if you want an all-NaN series with the same index:
                agg[f"{conc}_{st}_{t}"] = pd.Series(np.nan, index=df_vals.index)

# turn into DataFrame
agg_df = pd.DataFrame(agg)

# 4. row-wise Z-score -----------------------------------------------------
# (avoid zero-division by filling constant rows with zero)
agg_z = agg_df.sub(agg_df.mean(axis=1), axis=0) \
             .div(agg_df.std(axis=1).replace(0, np.nan), axis=0) \
             .fillna(0.0)

# 5. reorder columns for nicer display -----------------------------------
ordered_cols = [
    f"{c}_{s}_{t}"
    for c in concs
    for s in status
    for t in times_by_conc[c]
]
agg_z = agg_z.loc[:, [c for c in ordered_cols if c in agg_z.columns]]

# 6. heat-map -------------------------------------------------------------
sns.set_theme(context="paper", style="white")
vmax = np.nanmax(np.abs(agg_z.to_numpy()))

cg = sns.clustermap(
    agg_z,
    cmap="RdBu_r",
    vmin=-vmax, vmax=vmax,
    row_cluster=True,
    col_cluster=False,
    yticklabels=False,
    figsize=(6, 14),
    dendrogram_ratio=0.15,
    cbar_pos=(0.94, 0.25, 0.02, 0.5),
    cbar_kws=dict(label="Halo-zone intensity (Z-score)")
)

cg.ax_heatmap.set_xticklabels(
    cg.ax_heatmap.get_xmajorticklabels(), rotation=45, ha="right"
)
cg.ax_heatmap.set_xlabel("Concentration  •  Transformation status  •  Time-point")

cg.savefig("halo_heatmap_25mM_12.5mM.png", dpi=300, bbox_inches="tight")
plt.show()
