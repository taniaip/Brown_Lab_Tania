#!/usr/bin/env python3
"""
make_manhattan.py  -  draw a publication-style Manhattan plot

For Halo:
input : gwas_z24h.assoc.linear   (PLINK 1.9 output)
output: manhattan_z24h_pub_new_GWAS.png

For Clearing
input : gwas_z48h.assoc.linear   (PLINK 1.9 output)
output: manhattan_z48h_pub_new_GWAS.png
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# params
# assoc_file  = "gwas_z24h.assoc.linear"
# out_png     = "manhattan_z24h_pub_new_GWAS.png"
assoc_file  = "gwas_z48h.assoc.linear"
out_png     = "manhattan_clearing_pub_new_GWAS.png"
# p_cutoff    = 5e-19          # genome-wide line  #for halo
# mark_cutoff = 5e-19          # mark hits above this   #for halo
p_cutoff    = 1e-10           # genome-wide line  #for halo
mark_cutoff = 1e-10           # mark hits above this   #for halo
dot_size    = 7
top_dot_size = 30

# load + tidy
g = pd.read_csv(assoc_file, sep=r'\s+', usecols=["CHR", "BP", "SNP", "P"])

# map each unique CHR string to an integer 1..N
chr_map = {name: i+1 for i, name in enumerate(sorted(g["CHR"].unique()))}
g["CHR"] = g["CHR"].map(chr_map)

# cumulative position as before
chr_lengths = g.groupby("CHR")["BP"].max().sort_index()
offset = chr_lengths.cumsum().shift(fill_value=0)
g["pos_cum"] = g["BP"] + g["CHR"].map(offset)

# plot
plt.figure(figsize=(12, 5))
colors = ["#1f77b4", "#ffbf5e"]      # blue, tan

for i, (chrom, df) in enumerate(g.groupby("CHR")):
    plt.scatter(df["pos_cum"], -np.log10(df["P"]),
                s=dot_size, color=colors[i % 2])

# genome-wide significance line
plt.axhline(-np.log10(p_cutoff), color="red", linestyle="--", linewidth=1)

# highlight top hits
hits = g[g["P"] < mark_cutoff]
plt.scatter(hits["pos_cum"], -np.log10(hits["P"]),
            s=top_dot_size, color=colors[1], marker="D", edgecolor="k")

# chromosome labels centred in each block
xticks = (offset + chr_lengths / 2).values
plt.xticks(xticks, chr_lengths.index)
plt.xlabel("Chromosome")
plt.ylabel(r"$-\log_{10}(P)$")
# plt.title("Manhattan plot – z$_{24h}$ halo GWAS")
plt.title("Manhattan plot – z$_{48h}$ clearing GWAS")
plt.tight_layout()
plt.savefig(out_png, dpi=300)
print(f"✓ {out_png} written")
