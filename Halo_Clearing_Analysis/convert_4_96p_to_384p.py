#!/usr/bin/env python3
"""
Combine four 96-well blocks (B–M columns) into one 384-well plate.

Block positions in plates_updated.xlsx  (1-based Excel):
    Plate 1 : rows  5-12 , cols B-M
    Plate 2 : rows 17-24 , cols B-M
    Plate 3 : rows 29-36 , cols B-M
    Plate 4 : rows 41-48 , cols B-M
Mapping:
    P1 → odd rows, odd cols   (0,0)
    P2 → odd rows, even cols  (0,1)
    P3 → even rows, odd cols  (1,0)
    P4 → even rows, even cols (1,1)
"""
from pathlib import Path
import string
import pandas as pd

INPUT_FILE  = Path("plates_updated.xlsx")
OUTPUT_FILE = Path("combined_384_plate.xlsx")

ROW96   = list(string.ascii_uppercase[:8])      # A-H (8 rows)
COL96   = list(range(1, 13))                    # 1-12 (12 columns)
ROW384  = list(string.ascii_uppercase[:16])     # A-P
COL384  = list(range(1, 25))                    # 1-24

# (row0,row1,col0,col1) in **0-based iloc**
BLOCKS = [
    (4, 11, 1, 12),   # Plate 1  rows 5-12
    (16, 23, 1, 12),  # Plate 2  rows 17-24
    (28, 35, 1, 12),  # Plate 3  rows 29-36
    (40, 47, 1, 12),  # Plate 4  rows 41-48
]
OFFSETS = [(0,0), (0,1), (1,0), (1,1)]         # row_off, col_off

big = pd.read_excel(INPUT_FILE, sheet_name=0, header=None, engine="openpyxl")

plates96 = []
for (r0, r1, c0, c1) in BLOCKS:
    blk = big.iloc[r0:r1+1, c0:c1+1].copy()
    blk.index   = ROW96                     # A-H
    blk.columns = COL96                    # 1-12
    plates96.append(blk)

print("✓ extracted all four 96-well blocks")

# empty 384 DataFrame 
plate384 = pd.DataFrame(index=ROW384, columns=COL384, dtype=object)

# map each 96-well into the 384 lattice 
for (idx, (row_off, col_off)) in enumerate(OFFSETS):
    src = plates96[idx]
    for r96, row_lab in enumerate(ROW96):            # 0-7
        for c96 in COL96:                            # 1-12
            val = src.loc[row_lab, c96]
            if pd.isna(val):
                continue
            r384 = 2*r96 + row_off
            c384 = 2*(c96 - 1) + col_off            # subtract 1, then double
            plate384.iat[r384, c384] = val


with pd.ExcelWriter(OUTPUT_FILE, engine="openpyxl") as wrt:
    plate384.to_excel(wrt, sheet_name="384_plate", index_label="Row/Col")

print("✓ wrote combined 384-well plate →", OUTPUT_FILE)
