import pandas as pd
import os

# Load your expression data
df = pd.read_csv("final_data_annotated_merged_04052022.tab", encoding="latin1")

# Output directory (optional, will create if it doesn't exist)
outdir = "fusion_expression_files"
os.makedirs(outdir, exist_ok=True)

# For each gene, write a file with FID/IID/expression
genes = df['systematic_name'].unique()
for gene in genes:
    sub = df[df['systematic_name'] == gene][['Strain', 'tpm']].copy()
    sub['FID'] = sub['Strain']
    sub['IID'] = sub['Strain']
    sub = sub[['FID', 'IID', 'tpm']]
    # Clean gene name for file
    fname = os.path.join(outdir, f"{gene}_expression.txt")
    sub.to_csv(fname, sep='\t', index=False)

print(f"Done! Created {len(genes)} files in {outdir}/")
