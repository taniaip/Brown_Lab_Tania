import sys
import math
import pandas as pd
import matplotlib.pyplot as plt

# manually encode your DHY213-prlk120 baseline by (conc, time) in hours
# This is aquired form the results of this script:
# (base) tania@tanias-MacBook-Air analysis % python /Users/tania/Desktop/PETase\ project/Tania_imaging/2025-06-14/analysis/halo_data_reformat_original.py

# DHY213-120 baseline (per concentration and time):
# concentration  time
# 12.5           24h    -0.006710
#                48h    -0.011621
#                72h    -0.003487
#                8h     -0.019055
#                96h    -0.003144
# 25.0           24h     0.011836
#                48h    -0.066539
#                72h    -0.082399
#                8h     -0.026273
#                96h    -0.058981

# So DHY213-prlk120 baseline:
BASELINE = {
    (12.5, 8):  -0.019055,
    (12.5, 24): -0.006710,
    (12.5, 48): -0.011621,
    (12.5, 72): -0.003487,
    (12.5, 96): -0.003144,
    (25.0, 8):  -0.026273,
    (25.0, 24):  0.011836,
    (25.0, 48): -0.066539,
    (25.0, 72): -0.082399,
    (25.0, 96): -0.058981,
}

def plot_individual_and_atlas(filename):
    df = pd.read_excel(filename, engine='openpyxl')
    strains = df['strain'].tolist()
    times = [8, 24, 48, 72, 96]
    concs = [12.5, 25.0]
    
    # ---- 1) save each strain individually ----
    for strain in strains:
        row = df.loc[df['strain']==strain].iloc[0]
        plt.figure(figsize=(6,4))
        for conc, color in zip(concs, ('C0','C1')):
            xs, ys = [], []
            for t in times:
                cols = [f"pRLK139_{conc}_{t}h_{i}_colony.value" for i in (1,2,3)]
                raw = row[cols].astype(float).values
                norm = raw - BASELINE[(conc,t)]
                xs += [t]*3; ys += norm.tolist()
            plt.scatter(xs, ys, color=color, label=f"{conc} mM", alpha=0.7)
        plt.title(f"{strain}")
        plt.xticks(times); plt.xlabel("Time (h)"); plt.ylabel("Norm. intensity")
        plt.legend(fontsize='small'); plt.tight_layout()
        plt.savefig(f"{strain}_pRLK139_norm.png", dpi=150)
        plt.close()
    
    # ---- 2) build the atlas ----
    n = len(strains)
    cols = 7
    rows = math.ceil(n/cols)
    fig, axes = plt.subplots(rows, cols, figsize=(cols*2, rows*2), sharex=True, sharey=True)
    axes = axes.flatten()
    
    for ax, strain in zip(axes, strains):
        row = df.loc[df['strain']==strain].iloc[0]
        for conc, color in zip(concs, ('C0','C1')):
            xs, ys = [], []
            for t in times:
                cols_c = [f"pRLK139_{conc}_{t}h_{i}_colony.value" for i in (1,2,3)]
                raw = row[cols_c].astype(float).values
                norm = raw - BASELINE[(conc,t)]
                xs += [t]*3; ys += norm.tolist()
            ax.scatter(xs, ys, color=color, s=10, alpha=0.6)
        ax.set_title(strain, fontsize=6)
        ax.set_xticks(times); ax.tick_params(axis='both', labelsize=5)
    
    # turn off any unused subplots
    for ax in axes[n:]:
        ax.axis('off')
    
    fig.text(0.5, 0.04, 'Time (h)', ha='center')
    fig.text(0.04, 0.5, 'Norm. intensity', va='center', rotation='vertical')
    plt.tight_layout(rect=(0.05,0.05,1,1))
    fig.legend([plt.Line2D([],[],marker='o',linestyle=''),], ['12.5 mM / 25 mM'], 
               loc='upper right', title='Conc.', fontsize='small')
    atlas_name = "all_strains_prlk139_norm_atlas.png"
    fig.savefig(atlas_name, dpi=200)
    print(f"→ saved atlas {atlas_name}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python plot_prlk139_norm_atlas.py <excel_file>")
        sys.exit(1)
    plot_individual_and_atlas(sys.argv[1])
