import re

# User options
gff3_file = "Saccharomyces_cerevisiae.R64-1-1.112.gff3"
plink_prefix = "yeast1011"           # This should match your .bed/.bim/.fam prefix
window = 200                      # Window size 
output_script = "run_plink_windows.sh"

# Yeast chromosome mapping (roman to PLINK numeric, or as in .bim)
chr_map = {
    "I": "1", "II": "2", "III": "3", "IV": "4", "V": "5", "VI": "6", "VII": "7",
    "VIII": "8", "IX": "9", "X": "10", "XI": "11", "XII": "12", "XIII": "13",
    "XIV": "14", "XV": "15", "XVI": "16", "Mito": "17"  # Mitochondrial as 17 (change if needed)
}

plink_cmds = []

with open(gff3_file) as f:
    for line in f:
        if line.startswith("#"):
            continue
        cols = line.strip().split('\t')
        if len(cols) < 9:
            continue
        seqid, source, feature, start, end, score, strand, phase, attributes = cols
        if feature != "gene":
            continue
        # Only protein-coding genes
        if "biotype=protein_coding" not in attributes:
            continue
        # Extract gene name (ID=gene:GENENAME)
        match = re.search(r'ID=gene:([^;]+)', attributes)
        if not match:
            continue
        gene = match.group(1)
        chr_num = chr_map.get(seqid)
        if not chr_num:
            continue  # Skip non-standard chromosomes
        start = int(start)
        end = int(end)
        from_bp = max(1, start - window)
        to_bp = end + window
        cmd = (
            f"plink --bfile {plink_prefix} "
            f"--chr {chr_num} --from-bp {from_bp} --to-bp {to_bp} "
            f"--make-bed --out {gene}_window"
        )
        plink_cmds.append(cmd)

# Write all commands to a bash script
with open(output_script, "w") as f:
    f.write("#!/bin/bash\n\n")
    for cmd in plink_cmds:
        f.write(cmd + "\n")

print(f"Generated {len(plink_cmds)} PLINK commands in {output_script}")
