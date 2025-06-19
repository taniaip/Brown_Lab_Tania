#!/usr/bin/env python3
"""
PETase halo & clearing z-score finder

Outputs one Excel file with all z-scores and hit flags of ALL strains(unlike the list_all_hits.py 
that output just the ones taht are a hit at least in one timepoint/concnetration)
"""

from pathlib import Path
import re
import numpy as np
import pandas as pd

# CONFIG
INPUT_HALO   = Path("all_halo_merged_25mM.csv")        # 24 h + 42 h
INPUT_CLEAR  = Path("clearing_all_merged_12.5mM.csv")  # 48 h
OUTPUT_XLSX  = Path("zScore_all_strains_summary.xlsx")

Z_HIGH = 2.0   # z > 2 hit for halo plates
Z_LOW  = -2.0  # z < −2 hit for clearing plate
# helpers
def add_stats(df: pd.DataFrame, cols: list[str], label: str, tail: str):
    """
    Add mean_<label>, z_<label>, hit_<label> columns to df and return df slice.
    tail = "high" → hit if z > Z_HIGH
           "low"  → hit if z < Z_LOW
    """
    df = df.copy()
    df[cols] = df[cols].apply(pd.to_numeric, errors="coerce")

    mean_col = f"mean_{label}"
    df[mean_col] = df[cols].mean(axis=1, skipna=True)

    pop_mean = df[mean_col].mean(skipna=True)
    pop_std  = df[mean_col].std(ddof=0, skipna=True)
    z_col = f"z_{label}"
    df[z_col] = (df[mean_col] - pop_mean) / pop_std

    hit_col = f"hit_{label}"
    if tail == "high":
        df[hit_col] = df[z_col] > Z_HIGH
    else:
        df[hit_col] = df[z_col] < Z_LOW

    return df[["strain", mean_col, z_col, hit_col]]

# load data
halo   = pd.read_csv(INPUT_HALO)
clear  = pd.read_csv(INPUT_CLEAR)

# ---------- detect replicate columns -------------------------------------
pat_24h = re.compile(r"^X\d+\.colony\.value\.24h$")
pat_42h = re.compile(r"^X\d+\.colony\.value\.42h$")
pat_48h = re.compile(r"^X\d+\.colony\.value")    # all replicates

cols_24h = [c for c in halo.columns  if pat_24h.match(c)]
cols_42h = [c for c in halo.columns  if pat_42h.match(c)]
cols_48h = [c for c in clear.columns if pat_48h.match(c)]

# ---------- compute means, z-scores, hits ---------------------------------
sum_24 = add_stats(halo[["strain"]+cols_24h], cols_24h, "24h", "high")
sum_42 = add_stats(halo[["strain"]+cols_42h], cols_42h, "42h", "high")
sum_48 = add_stats(clear[["strain"]+cols_48h], cols_48h, "48h", "low")

# ---------- merge side-by-side -------------------------------------------
summary = (sum_24
           .merge(sum_42, on="strain", how="outer")
           .merge(sum_48, on="strain", how="outer"))

# ---------- keep only strains with ≥1 hit --------------------------------
hit_any = (
    summary["hit_24h"].fillna(False)
    | summary["hit_42h"].fillna(False)
    | summary["hit_48h"].fillna(False)
)
summary = summary.sort_values("strain").reset_index(drop=True)

# ---------- write Excel ---------------------------------------------------
with pd.ExcelWriter(OUTPUT_XLSX, engine="openpyxl") as wrt:
    summary.to_excel(wrt, sheet_name="hits", index=False)

