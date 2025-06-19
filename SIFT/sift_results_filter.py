import pandas as pd

# File name
file_name = "SIFT4G_results_clearing_over_9/clearing_hits_over_9_SIFTannotations.xlsx"
output_sheet = "Deleterious"

# Read the data (works for .xls or .xlsx; if issues, save as .xlsx and change file extension)
df = pd.read_excel(file_name)

# Filter for DELETERIOUS predictions
deleterious_df = df[
    (df['SIFT_PREDICTION'].str.upper() == 'DELETERIOUS') |
    (df['SIFT_PREDICTION'].str.upper() == 'DELETERIOUS (*WARNING! LOW CONFIDENCE)')
]

# Write the filtered data to a new sheet in the same file
with pd.ExcelWriter(file_name, mode='a', engine='openpyxl', if_sheet_exists='replace') as writer:
    deleterious_df.to_excel(writer, sheet_name=output_sheet, index=False)

print("✓ Filtered sheet 'Deleterious' added to", file_name)
