#!/usr/bin/env python3
"""
Identify genes whose natural loss-of-function (LOF) mutations are
significantly enriched among high-halo (Z≥1.5) strains.
Output:  LOF_enriched_genes_highZ.xlsx
"""

import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# 1. Load LOF matrix
lof_raw = pd.read_excel("1011_LOF_Peter_2018.xlsx", engine="openpyxl")
gene_col = "NAME"  # Ensure this matches your file's gene name column

# Set gene name as index; select only strain columns (skip metadata)
strain_columns = lof_raw.columns[7:]  # all columns after STOP
lof_bin = lof_raw.set_index(gene_col)[strain_columns]
lof_bin = (lof_bin == 1).astype(int)  # 1 = LOF, 0 = WT

# Uppercase for robustness
lof_bin.index = lof_bin.index.astype(str).str.upper().str.strip()
lof_bin.columns = lof_bin.columns.astype(str).str.strip()

# 2. Load high-halo strains
z_df = pd.read_excel("zScore_hits_summary.xlsx")
z_df["strain"] = z_df["strain"].str.strip()
# Identify strains that are 'high' in 24h or 42h
high = set(z_df.loc[(z_df["hit_24h"]==True) | (z_df["hit_42h"]==True), "strain"])
all_set = set(lof_bin.columns)
high &= all_set  # intersection
non_high = list(all_set - high)
high_list = list(high)

# 3. Fisher exact test per gene
records = []
for g in lof_bin.index:
    n11 = lof_bin.loc[g, high_list].sum()      # LOF & high-halo
    n10 = lof_bin.loc[g, non_high].sum()       # LOF & non-high
    n01 = len(high_list) - n11                 # WT  & high-halo
    n00 = len(non_high)  - n10                 # WT  & non-high
    odds, p = fisher_exact([[n11, n01], [n10, n00]], alternative="greater")
    records.append((g, n11, n10, odds, p))

res_df = pd.DataFrame(
    records, columns=["Gene", "LOF_high", "LOF_rest", "OddsRatio", "P"]
)

# 4. Multiple-test correction
res_df["FDR"] = multipletests(res_df["P"], method="fdr_bh")[1]
sig = (res_df[(res_df["FDR"] < 0.05) & (res_df["OddsRatio"] > 1)]
         .sort_values("OddsRatio", ascending=False))

# 5. Export
out_file = "LOF_enriched_genes_highZ.xlsx"
sig.to_excel(out_file, index=False, engine="openpyxl")
print(f"{len(sig)} genes significantly enriched (FDR < 0.05).")
print(f"Results written to {out_file}")
