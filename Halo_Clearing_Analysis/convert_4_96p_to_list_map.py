#!/usr/bin/env python3
"""
convert_4_96p_to_list_map.py
Read `combined_384_plate.xlsx`  ►  write `plate_map.xlsx`

Output columns
p   - plate number (always 1 here)
c   - column  (1 … 24)
r   - row     (1 … 16, where A→1 … P→16)
gene - strain / gene code from the cell
"""

import pandas as pd # type: ignore
from pathlib import Path
import string

# paths 
INPUT_XLSX  = Path("series2_plates_updated.xlsx")     # source file
OUTPUT_XLSX = Path("plate_map_series2.xlsx")         # result

# load the plate 
# if the sheet is named '384_plate'; else sheet_name=0 reads the first sheet
plate = pd.read_excel(INPUT_XLSX, sheet_name=0, index_col=0)

# ensure we have exactly 16 rows (A-P) and 24 columns (1-24)
row_lookup = {letter: i+1 for i, letter in enumerate(string.ascii_uppercase[:16])}

records = []
for row_label, row in plate.iterrows():
    r = row_lookup.get(str(row_label).strip().upper())
    if r is None:
        continue                           # skip unexpected row labels
    for c, gene in row.items():
        if pd.isna(gene):
            continue                       # skip empty wells
        records.append({"p": 1, "c": c, "r": r, "gene": str(gene)})

tidy = pd.DataFrame(records).sort_values(["r", "c"]).reset_index(drop=True)

# write 
tidy.to_excel(OUTPUT_XLSX, index=False)
print(f"✓ wrote {len(tidy)} rows to {OUTPUT_XLSX}")
