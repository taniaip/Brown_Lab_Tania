#!/usr/bin/env python3
"""
growth_stats.py
Analyse 768-well PETase growth metrics.

Input  : 768_growth_full_metrics.xlsx   (first sheet)
Output : significant_growth_effects.xlsx
         ├─ Sheet “120_vs_139”   - strains whose 120 vs 139 t_gen differ (P<0.05)
         └─ Sheet “vs_DHY213”    - strains that differ from DHY213 with the same plasmid (P<0.1)

"""

import re, sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

# CONFIG
IN_FILE   = Path("768_growth_full_metrics.xlsx")
OUT_FILE  = Path("significant_growth_effects.xlsx")
SIGMA_TAIL = 0.05        # discard wells in the lowest 5 % of σ
P_CUTOFF   = 0.1         # significance threshold for t-tests
MIN_N      = 2           # require at least this many wells per group
# helpers
def parse_label(label: str):
    """
    Split 'APV 120-2'  →  base='APV', plasmid='120'
           'DHY213 - 139' → base='DHY213', plasmid='139'
    """
    m = re.search(r"(120|139)", label)
    if not m:
        return label.strip(), None
    plasmid = m.group(1)
    base = label[:m.start()].strip(" -_")
    base = re.sub(r"[-\s]\d+$", "", base)   # drop trailing '-1', ' 3' etc.
    return base, plasmid
# load & clean
df = pd.read_excel(IN_FILE, sheet_name=0)

required = {"strain", "t_gen", "sigma"}
if not required.issubset(df.columns):
    sys.exit(f"ERROR: {IN_FILE} is missing columns {required - set(df.columns)}")

# add 'base' and 'plasmid' columns
df[["base", "plasmid"]] = df["strain"].apply(
    lambda s: pd.Series(parse_label(str(s)))
)

# drop rows without a recognised plasmid
df = df.dropna(subset=["plasmid"])

# 1) throw out wells whose σ is “too small”
sigma_mean = df["sigma"].mean(skipna=True)
sigma_std  = df["sigma"].std(skipna=True)
sigma_cut  = sigma_mean - abs(np.percentile(df["sigma"].dropna(), SIGMA_TAIL*100))
# equivalently: sigma_mean - z_0.95 * std  (lower 5 %)
clean = df[df["sigma"] >= sigma_cut].copy()

# function to run Welch t-test
def welch(group_a, group_b):
    if len(group_a) < MIN_N or len(group_b) < MIN_N:
        return np.nan, np.nan, np.nan
    t, p = ttest_ind(group_a, group_b, equal_var=False, nan_policy="omit")
    return p, group_a.mean(), group_b.mean()
# analysis 1: 120 vs 139
rows_120_139 = []
for base, grp in clean.groupby("base"):
    tg120 = grp.loc[grp["plasmid"]=="120", "t_gen"]
    tg139 = grp.loc[grp["plasmid"]=="139", "t_gen"]
    p, mean120, mean139 = welch(tg120, tg139)
    if np.isfinite(p) and p < P_CUTOFF:
        rows_120_139.append({
            "strain"         : base,
            "n_120"          : len(tg120),
            "n_139"          : len(tg139),
            "mean_t_gen_120" : mean120,
            "mean_t_gen_139" : mean139,
            "ratio_120/139"  : mean120 / mean139 if mean139 else np.nan,
            "p_value"        : p,
        })
sheet1 = pd.DataFrame(rows_120_139).sort_values("ratio_120/139")

# analysis 2: each strain vs DHY213 (same plasmid)
rows_vs = []
for plasmid, grp_pl in clean.groupby("plasmid"):
    dhy = grp_pl.loc[grp_pl["base"]=="DHY213", "t_gen"]
    for base, grp in grp_pl.groupby("base"):
        if base == "DHY213":
            continue
        tg   = grp["t_gen"]
        p, mean_dhy, mean_strain = welch(dhy, tg)
        if np.isfinite(p) and p < P_CUTOFF:
            rows_vs.append({
                "strain"            : base,
                "plasmid"           : plasmid,
                "n_strain"          : len(tg),
                "n_DHY213"          : len(dhy),
                "mean_t_gen_strain" : mean_strain,
                "mean_t_gen_DHY213" : mean_dhy,
                "ratio_DHY/strain"  : mean_dhy / mean_strain if mean_strain else np.nan,
                "p_value"           : p,
            })
sheet2 = pd.DataFrame(rows_vs).sort_values(["plasmid","ratio_DHY/strain"])

# write workbook
with pd.ExcelWriter(OUT_FILE, engine="openpyxl") as xls:
    sheet1.to_excel(xls, sheet_name="120_vs_139", index=False)
    sheet2.to_excel(xls, sheet_name="vs_DHY213", index=False)

print("✓ significant_growth_effects.xlsx created")
