import pandas as pd
from statsmodels.stats.multitest import multipletests

# 1) Read the merged results
# Adjust sep="\t" if your .dat is tab-delimited, or sep="\s+" for whitespace
df = pd.read_csv("FUSION_ALL_CHR.dat", sep="\t", dtype=str)

# 2) Convert the P column to numeric
df["TWAS.P"] = pd.to_numeric(df["TWAS.P"], errors="coerce")

# 3) Compute FDR (Benjamini-Hochberg)
# Fill NaNs with 1.0 so multipletests returns an array of the same length
pvals = df["TWAS.P"].fillna(1.0)
_, fdr, _, _ = multipletests(pvals, method="fdr_bh")
df["FDR"] = fdr

# 4) Identify genes passing FDR < 0.05
sig = df[df["FDR"] < 0.05]
print(f"Number of significant genes: {len(sig)}")

# 5) Write them out to Excel
output_file = "FUSION_significant_genes.xlsx"
with pd.ExcelWriter(output_file, engine="openpyxl") as writer:
    sig.to_excel(writer, sheet_name="Significant_Genes", index=False)

print(f"Significant genes written to {output_file}")

