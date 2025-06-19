#!/usr/bin/env python3
#python txt_to_vcf.py --input corrected_halo_genomewide_hits_pos_beta.txt --fasta Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa --output halo_hits.vcf
"""
Convert GWAS text hits (with columns CHR, BP, A1, ...) to a minimal VCF file.
Usage:
    python text_to_vcf.py --input halo_pos_beta_hits.txt --fasta reference.fa --output hits.vcf
"""

import argparse
import pandas as pd
import re
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(description="Convert GWAS hits (text) to minimal VCF.")
    parser.add_argument('--input', '-i', required=True, help='Input hits text file (whitespace-delimited)')
    parser.add_argument('--fasta', '-f', required=True, help='Reference FASTA file')
    parser.add_argument('--output', '-o', required=True, help='Output VCF file')
    args = parser.parse_args()

    # 1. Load hits table
    hits = pd.read_csv(args.input, delim_whitespace=True)

    # 2. Load reference FASTA (keys: I, II, ...)
    fasta = {record.id: str(record.seq).upper() for record in SeqIO.parse(args.fasta, "fasta")}

    # 3. Map chromosomeN → Roman numeral
    roman = {str(i): r for i, r in enumerate(
        ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"], 1)}

    # 4. Output VCF
    with open(args.output, 'w') as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for _, row in hits.iterrows():
            chrstr = str(row['CHR'])
            m = re.match(r'chromosome(\d+)', chrstr)
            if not m:
                continue
            chrom = roman.get(m.group(1))
            pos = int(row['BP'])
            alt = str(row['A1'])
            # Check and get REF base
            if chrom not in fasta or pos-1 >= len(fasta[chrom]):
                continue
            ref = fasta[chrom][pos-1]
            snp_id = str(row['SNP']) if 'SNP' in hits.columns else '.'
            # Only keep SNPs (not indels)
            if len(ref) == 1 and len(alt) == 1:
                f.write(f"{chrom}\t{pos}\t{snp_id}\t{ref}\t{alt}\t.\t.\t.\n")

    print(f"✓ Wrote VCF: {args.output}")

if __name__ == '__main__':
    main()
