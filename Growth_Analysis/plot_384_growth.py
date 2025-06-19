#!/usr/bin/env python3
"""
plot_384_growth.py
Creates a single-page PDF with all 384 growth curves + strain labels.

Inputs--
• 15_05_25_Tania_validation_ps1_A_384_T_23.txt   - Tecan export
• combined_384_plate.xlsx  (sheet '384_plate')    - strain lookup

Outputs---
• growth_curves_384.pdf   - 16x24 grid of fitted curves
• growth_metrics.csv      - k, n0, r, t_mid, doubling_h, sigma …
"""
import io, sys
from pathlib import Path
import numpy as np
import pandas as pd # type: ignore
from scipy.optimize import curve_fit # type: ignore
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# RAW_TXT  = Path("15_05_25_Tania_validation_ps1_A_384_T_23.txt")
RAW_TXT  = Path("15_05_25_tania_validation#2_ps1_A_384_T_24.txt")

LAYOUT   = Path("combined_384_plate.xlsx")
OUT_PDF  = Path("growth_curves_384.pdf")
OUT_CSV  = Path("growth_metrics.csv")

READ_INT_H = 0.5            # 30-min reads
ROLL_MED   = 5              # smoothing for plotting raw OD
MAXFEV     = 10000          # optimiser iterations
#-

def parse_tecan(path: Path) -> pd.DataFrame:
    """Return dataframe of first [Data:] block."""
    lines = path.read_text().splitlines()
    blk, keep = [], False
    for ln in lines:
        if ln.startswith("[Data:]"):
            keep = True
            continue
        if keep:
            if ln.startswith("[") and not ln.startswith("[Data:]"):
                break
            blk.append(ln)
    if blk[0].startswith("Read = "):
        blk[0] = blk[0].replace("Read = ", "", 1)
    return pd.read_csv(io.StringIO("\n".join(blk)), sep="\t")

def logistic(t, n0, k, r, t_mid):
    return n0 + (k - n0) / (1 + np.exp(-r * (t - t_mid)))

def fit_curve(t, y):
    """Fit 4-param logistic; return k,n0,r,t_mid,sigma."""
    y = np.clip(y, 1e-3, None)
    lb = [0, 0,   0,   0]
    ub = [3, 3,  10,  t[-1]]
    p0 = [y.min(), y.max(), 0.5, t[len(t)//2]]
    try:
        popt, _ = curve_fit(logistic, t, y, p0=p0,
                            bounds=(lb, ub), maxfev=MAXFEV)
        resid = y - logistic(t, *popt)
        sigma = np.sqrt(np.mean(resid**2))
        return *popt, sigma
    except Exception:
        return (np.nan,)*5

def main():
    # load raw OD
    raw = parse_tecan(RAW_TXT)
    time_h = (raw["RunT[s]"].to_numpy()/3600 if "RunT[s]" in raw.columns
              else np.arange(len(raw))*READ_INT_H)
    well_cols = [c for c in raw.columns if c.startswith("Well_")]

    # load strain names -
    layout = (pd.read_excel(LAYOUT, sheet_name="384_plate", index_col=0)
              .stack().reset_index().rename(columns={0:"strain"}))
    layout["well"] = (layout["Row/Col"]
                      + layout["level_1"].astype(int).map("{:02d}".format))
    
    strain_dict = dict(zip(layout["well"], layout["strain"]))
    row_major = [c for c in raw.columns if c.startswith("Well_")]
    well_cols = [row_major[row + col*16]        # row 0-15, col 0-23
                for col in range(24)
                for row in range(16)]
    # prepare PDF canvas 
    pdf = PdfPages(OUT_PDF)
    fig_w, fig_h = 24, 16   # inches; 1×1 per mini-plot
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = fig.add_gridspec(16, 24, wspace=0.05, hspace=0.05)

    metrics = []
    for idx, col in enumerate(well_cols, start=1):
        row_i = (idx-1)//24
        col_i = (idx-1)%24
        ax = fig.add_subplot(gs[row_i, col_i])

        y_raw = raw[col].to_numpy()
        k,n0,r,t_mid,sigma = fit_curve(time_h, y_raw)
        doubling = np.log(2)/r if np.isfinite(r) and r>0 else np.nan
        well_id = "ABCDEFGHIJKLMNOP"[row_i] + f"{col_i+1:02d}"
        strain  = strain_dict.get(well_id, "")

        # plot raw (smoothed) & fit
        y_smooth = (pd.Series(y_raw)
                    .rolling(ROLL_MED, center=True, min_periods=1)
                    .median())
        ax.plot(time_h, y_smooth, '.', ms=2, color="grey")
        # if np.isfinite(k):
        #     ax.plot(time_h, logistic(time_h, n0,k,r,t_mid),
        #             lw=0.8, color="#1f77b4")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlim(0, time_h[-1]); ax.set_ylim(0, max(1.2, y_raw.max()*1.1))
        ax.text(0.02, 0.95, well_id, transform=ax.transAxes,
                fontsize=6, va='top', ha='left')
        ax.text(0.98, 0.95, strain, transform=ax.transAxes,
                fontsize=5, va='top', ha='right')

        metrics.append(dict(sample=well_id, strain=strain,
                            k=k, n0=n0, r=r, t_mid=t_mid,
                            doubling_h=doubling, sigma=sigma))

    fig.suptitle("Growth curves - 384-well plate", fontsize=14)
    pdf.savefig(fig, dpi=300, bbox_inches="tight")
    pdf.close()
    plt.close(fig)
    print(f"✓ {OUT_PDF.name} written (all 384 curves on one page)")

    # save metrics table 
    pd.DataFrame(metrics).to_csv(OUT_CSV, index=False)
    print(f"✓ {OUT_CSV.name} written")

if __name__ == "__main__":
    try: main()
    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)
