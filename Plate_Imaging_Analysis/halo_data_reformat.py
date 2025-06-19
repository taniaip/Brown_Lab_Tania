#!/usr/bin/env python3
# halo_reformat_minimal.py  ─  stripped-down version
# -------------------------------------------------------------
#  input : halo.data.med.csv
#  output: halo_clearing_cleaned_data.xlsx
# -------------------------------------------------------------
import pandas as pd
import re
pd.options.mode.chained_assignment = "warn"   # keep warnings, not errors

# 1. load once 
df_all = pd.read_csv("halo.data.med.csv")

# 2. parse metadata ON THE FULL TABLE  (key fix)─
df_all["base_strain"]  = df_all["strain"].str.replace(r"-\d+$", "", regex=True)
df_all["plasmid"]      = df_all["strain"].str.endswith("139") \
                           .map({True: "pRLK139", False: "pRLK120"})
df_all["concentration"] = df_all["image.name"].str.extract(r"(\d+(?:\.\d+)?)mM")[0].astype(float)
df_all["time"]          = df_all["image.name"].str.extract(r"_(\d+h)")[0]

# 3. split DHY213 baseline from the rest 
base = df_all[
    (df_all["strain"].str.startswith("DHY213")) &
    (df_all["plasmid"] == "pRLK120")            # DHY213 with the empty vector
]
df   = df_all[~df_all["strain"].str.startswith("DHY213")].copy()

# 4. replicate index 1…n in each block 
df = df.sort_values(["base_strain", "plasmid", "concentration", "time"])
df["rep"] = (df.groupby(["base_strain", "plasmid", "concentration", "time"])
               .cumcount() + 1)

# tidy label e.g. pRLK139_25.0_48h_2_colony.value
df["label"] = (df["plasmid"] +"_"+ df["concentration"].astype(str) +"_"+
               df["time"] +"_"+ df["rep"].astype(str)+"_colony.value")

# 5. pivot replicates 
rep_wide = df.pivot(index="base_strain",
                    columns="label",
                    values="colony.value")

# --- 6. block mean & DHY213-normalised mean including timepoint 
ave_raw = (df.groupby(["base_strain", "plasmid", "concentration", "time"])
             ["colony.value"]
             .mean()
             .unstack(["plasmid", "concentration", "time"]))

# DHY213 baseline: pRLK120 ONLY, one value per concentration *and* time
baseline = (base.groupby(["concentration", "time"])["colony.value"]
              .mean())

print("\nDHY213-120 baseline (per concentration and time):")
print(baseline)

# subtract that baseline from every block that shares the same conc & time
norm_raw = ave_raw.copy()
for (conc, t), b in baseline.items():  # conc = 12.5, 25.0; t = 24h, 96h, etc.
    norm_raw.loc[:, pd.IndexSlice[:, conc, t]] = (
        ave_raw.loc[:, pd.IndexSlice[:, conc, t]] - b
    )

# 3. rename columns after the subtraction, now with time
ave_int       = ave_raw.copy()
ave_int.columns = [f"{p}_{c}_{t}_avg"       for p, c, t in ave_raw.columns]

norm_ave_int  = norm_raw.copy()
norm_ave_int.columns = [f"{t}_{p}_{c}_avg_norm" for p, c, t in norm_raw.columns]

# 7. save 
out = pd.concat([rep_wide, ave_int, norm_ave_int], axis=1)
out.index.name = "strain"
out.reset_index().to_excel("halo_clearing_cleaned_data.xlsx", index=False)

print("✓ halo_clearing_cleaned_data.xlsx written – DHY213 rows excluded")
