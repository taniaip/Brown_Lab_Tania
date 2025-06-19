#!/usr/bin/env python3
"""
make_full_metrics.py
Produces: sample  k  n0  r  t_mid  t_gen  auc_l  auc_e  sigma  note  strain
"""

import io, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.integrate import trapezoid

RAW_TXT   = Path("15_05_25_tania_validation#2_ps1_A_384_T_24.txt")
# RAW_TXT   = Path("15_05_25_Tania_validation_ps1_A_384_T_23.txt")
LAYOUT_XL = Path("combined_384_plate.xlsx")          # sheet '384_plate'
OUT_XLSX  = Path("growth_full_metrics.xlsx")

READ_INT_H = 0.5          # 30-min cycle
MAX_ITER   = 10000

def parse_tecan_block(path: Path) -> pd.DataFrame:
    text = path.read_text().splitlines()
    in_block, blk = False, []
    for ln in text:
        if ln.startswith("[Data:]"):
            in_block = True
            continue
        if in_block:
            if ln.startswith("[") and not ln.startswith("[Data:]"):
                break
            blk.append(ln)
    if blk[0].startswith("Read = "):
        blk[0] = blk[0].replace("Read = ", "", 1)
    return pd.read_csv(io.StringIO("\n".join(blk)), sep="\t")

def logistic(t, n0, k, r, t_mid):
    return n0 + (k - n0) / (1.0 + np.exp(-r * (t - t_mid)))

def fit_one(t, y):
    """Return k,n0,r,t_mid,sigma (np.nan on fail)."""
    # basic bounds
    lb = [0, 0,     0,   0]
    ub = [3, 3,   10,  t[-1]]
    p0 = [y.min(), y.max(), 0.5, t[len(t)//2]]
    try:
        popt, pcov = curve_fit(logistic, t, y, p0=p0,
                               bounds=(lb, ub), maxfev=MAX_ITER)
        residuals = y - logistic(t, *popt)
        sigma = np.sqrt(np.mean(residuals**2))
        return *popt, sigma
    except Exception:
        return (np.nan, np.nan, np.nan, np.nan, np.nan)
# or, equivalently:
# return tuple([np.nan]*5)

def main():
    df_raw = parse_tecan_block(RAW_TXT)
    t = (df_raw["RunT[s]"].to_numpy()/3600.0
         if "RunT[s]" in df_raw.columns
         else np.arange(len(df_raw))*READ_INT_H)

    wells = [c for c in df_raw.columns if c.startswith("Well_")]
    rows = []
    for col in wells:
        y = df_raw[col].clip(lower=1e-3).to_numpy()
        if y.max() < 0.02:
            note = "flat"
            k = n0 = r = t_mid = t_gen = auc_l = auc_e = sigma = np.nan
        else:
            k, n0, r, t_mid, sigma = fit_one(t, y)
            t_gen = np.log(2)/r if r and not np.isnan(r) else np.nan
            auc_l = trapezoid(y, t)
            auc_e = trapezoid(np.log(y), t)
            note = ""
        idx = int(col.split("_")[1])
        rc  = "ABCDEFGHIJKLMNOP"[(idx-1)//24] + f"{(idx-1)%24+1:02d}"
        rows.append(dict(sample=rc, k=k, n0=n0, r=r, t_mid=t_mid,
                         t_gen=t_gen, auc_l=auc_l, auc_e=auc_e,
                         sigma=sigma, note=note))

    metrics = pd.DataFrame(rows)

    # merge strain names 
    layout = (pd.read_excel(LAYOUT_XL, sheet_name="384_plate", index_col=0)
              .stack().rename("strain").reset_index())
    layout["sample"] = (layout["Row/Col"] +
                        layout["level_1"].astype(int).map("{:02d}".format))
    metrics = metrics.merge(layout[["sample", "strain"]], on="sample", how="left")

    # reorder columns
    metrics = metrics[["sample","k","n0","r","t_mid","t_gen",
                       "auc_l","auc_e","sigma","note","strain"]]

    metrics.to_excel(OUT_XLSX, index=False)
    print("✓ wrote", OUT_XLSX)

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)
