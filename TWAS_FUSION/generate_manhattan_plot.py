
#!/usr/bin/env python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import csv

# 1) Try to sniff, but fallback to tab if that fails
with open("FUSION_ALL_CHR.dat", "r") as f:
    sample = f.read(2048)
try:
    dialect = csv.Sniffer().sniff(sample, delimiters="\t ")
    sep = dialect.delimiter
    print("Detected delimiter:", repr(sep))
except csv.Error:
    sep = "\t"
    print("Could not auto-detect delimiter, falling back to tab.")

# 2) Read everything as strings, strip whitespace from column names
df = pd.read_csv("FUSION_ALL_CHR.dat", sep=sep, dtype=str, engine="python")
df.columns = df.columns.str.strip()

# 3) Drop any column that is entirely missing or literally "NA"
#    (first replace any "NA" strings with actual np.nan, then drop all-nan columns)
df.replace("NA", np.nan, inplace=True)
df.dropna(axis=1, how="all", inplace=True)
print("Columns after dropping all-NA:", df.columns.tolist())

# 4) Identify the TWAS P‐value column (must start with "TWAS.P")
pcols = [c for c in df.columns if c.upper().startswith("TWAS.P")]
if not pcols:
    raise Exception("No column starting with 'TWAS.P' found in: " + ", ".join(df.columns))
pcol = pcols[0]
print("Using P-value column:", pcol)

# 5) Coerce the P‐value column to numeric, drop invalid rows
df[pcol] = pd.to_numeric(df[pcol], errors="coerce")
before = len(df)
df = df.dropna(subset=[pcol])
print("Kept %d / %d rows with valid P-values" % (len(df), before))

# 6) Compute –log10(P) and a simple position index
df["minus_log10_P"] = -np.log10(df[pcol].astype(float))
df["pos"] = np.arange(len(df))

# 7) Plot Manhattan
plt.figure(figsize=(10, 4))
for chr_val in sorted(df["CHR"].astype(str).unique()):
    sub = df[df["CHR"].astype(str) == chr_val]
    plt.scatter(sub["pos"], sub["minus_log10_P"], s=10, alpha=0.6, label=chr_val)

plt.xlabel("Gene index")
plt.ylabel("-log10(P-value)")
plt.title("TWAS Manhattan Plot")
plt.legend(title="Chr", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize="small")
plt.tight_layout()
plt.savefig("TWAS_manhattan.png", dpi=300)
print("Saved TWAS_manhattan.png")


