import argparse
import gffutils
import pprint
import pandas as pd
from Bio import SeqIO
from cyvcf2 import VCF

# Make printouts wider for debugging
pprint.pprint = lambda obj, **kwargs: pprint.PrettyPrinter(width=120).pprint(obj)

parser = argparse.ArgumentParser()
parser.add_argument("-g", "--genes", dest="genes_gff3_file", type=str, required=True, help="Input GFF3 file with gene information.")
parser.add_argument("-f", "--fasta", dest="fasta_file", type=str, required=True, help="Input FASTA file containing genomic sequences.")
parser.add_argument("-i", "--hits", dest="hits_file", type=str, required=True, help="Input SNP GWAS hits text file.")
parser.add_argument("-s", "--summary", dest="strain_summary_file", type=str, required=True, help="Input summary Excel file for strains.")
parser.add_argument("-v", "--vcf", dest="vcf_file", type=str, required=True, help="Input VCF file for strain variants.")
parser.add_argument("-o", "--output", dest="output_file", type=str, required=True, help="Output Excel file with annotations and strain matrix.")
args = parser.parse_args()

def load_fasta_sequences(fasta_file):
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

def reverse_complement(seq):
    base_convertor = {'T':'A', 'A':'T', 'C':'G', 'G':'C'}
    return ''.join(base_convertor.get(b, b) for b in seq[::-1])

codon_to_AA = {
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G'
}

english_to_roman = {
    "chromosome1": "I", "chromosome2": "II", "chromosome3": "III", "chromosome4": "IV", "chromosome5": "V",
    "chromosome6": "VI", "chromosome7": "VII", "chromosome8": "VIII", "chromosome9": "IX", "chromosome10": "X",
    "chromosome11": "XI", "chromosome12": "XII", "chromosome13": "XIII", "chromosome14": "XIV", "chromosome15": "XV", "chromosome16": "XVI"
}

# --- Annotation Section ---

# Load data
fasta = load_fasta_sequences(args.fasta_file)
db = gffutils.create_db(args.genes_gff3_file, dbfn=':memory:', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
hits = pd.read_csv(args.hits_file, sep=r'\s+', comment="#")
hits['gene_found'] = False
hits['gene_id'] = None
hits['strand'] = None
hits['gene_function'] = None
hits['mutation_type'] = None
hits['original_codon'] = None
hits['new_codon'] = None
hits['original_AA'] = None
hits['new_AA'] = None

for idx, row in hits.iterrows():
    chrom = row['CHR']
    BP = int(row['BP'])
    allele = row['A1']
    roman_chr = english_to_roman.get(chrom)
    if not roman_chr or roman_chr not in fasta:
        continue
    found = False
    for feature in db.region(seqid=roman_chr, start=BP, end=BP):
        if feature.featuretype == 'gene':
            hits.at[idx, 'gene_found'] = True
            hits.at[idx, 'gene_id'] = feature.attributes.get('ID', [''])[0]
            hits.at[idx, 'gene_function'] = feature.attributes.get('description', [''])[0] if 'description' in feature.attributes else ""
        if feature.featuretype == 'CDS' and feature.start <= BP <= feature.end:
            found = True
            strand = feature.strand
            cds_start = feature.start
            cds_end = feature.end
            codon_pos = (BP - cds_start) if strand == '+' else (cds_end - BP)
            codon_frame = codon_pos % 3
            if strand == '+':
                codon_start = BP - codon_frame
                codon_seq = fasta[roman_chr][codon_start-1:codon_start+2].upper()
                hits.at[idx, 'original_codon'] = codon_seq
                hits.at[idx, 'strand'] = strand
                hits.at[idx, 'original_AA'] = codon_to_AA.get(codon_seq, 'X')
                if len(allele) == 1:
                    codon_list = list(codon_seq)
                    codon_list[codon_frame] = allele
                    new_codon = ''.join(codon_list)
                    hits.at[idx, 'new_codon'] = new_codon
                    new_AA = codon_to_AA.get(new_codon, 'X')
                    hits.at[idx, 'mutation_type'] = 'point'
                    hits.at[idx, 'new_AA'] = new_AA
                else:
                    hits.at[idx, 'mutation_type'] = 'frameshift'
                    hits.at[idx, 'new_codon'] = None

            elif strand == '-':
                codon_start = BP + codon_frame
                codon_seq = fasta[roman_chr][codon_start-3:codon_start].upper()
                hits.at[idx, 'strand'] = strand
                orig_codon_coding = reverse_complement(codon_seq)
                hits.at[idx, 'original_codon'] = orig_codon_coding
                hits.at[idx, 'original_AA'] = codon_to_AA.get(orig_codon_coding, 'X')
                if len(allele) == 1:
                    codon_list = list(codon_seq)
                    codon_frame_corrected = codon_frame * -1 -1
                    codon_list[codon_frame_corrected] = allele
                    new_codon_plus = ''.join(codon_list)
                    new_codon_coding = reverse_complement(new_codon_plus)
                    hits.at[idx, 'new_codon'] = new_codon_coding
                    new_AA = codon_to_AA.get(new_codon_coding, 'X')
                    hits.at[idx, 'mutation_type'] = 'point'
                    hits.at[idx, 'new_AA'] = new_AA
                else:
                    hits.at[idx, 'mutation_type'] = 'frameshift'
                    hits.at[idx, 'new_codon'] = None
            break
    if not hits.at[idx, 'gene_found']:
        hits.at[idx, 'mutation_type'] = 'intergenic'

# --- Presence/Absence Matrix Section ---

strain_df = pd.read_excel(args.strain_summary_file)
strain_list = list(strain_df['strain'].unique())
vcf = VCF(args.vcf_file)
vcf_strains = list(vcf.samples)
strain_list_in_vcf = [s for s in strain_list if s in vcf_strains]
if not strain_list_in_vcf:
    raise ValueError("No overlapping strains between summary and VCF!")

# Prepare columns for each strain (init with 0)
for strain in strain_list_in_vcf:
    hits[strain] = 0

# Prepare a lookup: For each hit, get presence in each strain
snp_lookup = {}
for idx, row in hits.iterrows():
    snp_id = f"{row['CHR']}_{int(row['BP'])}_{row['A1']}"
    if snp_id not in snp_lookup:
        snp_lookup[snp_id] = []
    snp_lookup[snp_id].append(idx)

for variant in vcf:
    for alt in variant.ALT:
        snp_id = f"{variant.CHROM}_{variant.POS}_{alt}"
        if snp_id in snp_lookup:
            for idx in snp_lookup[snp_id]:
                for strain in strain_list_in_vcf:
                    s_idx = vcf_strains.index(strain)
                    gt = variant.genotypes[s_idx][:2]
                    present = int(any(a == 1 for a in gt if a is not None))
                    hits.at[idx, strain] = present

# Save one file with both annotation and matrix
hits.to_excel(args.output_file, index=False)
print(f"✓ Done! Output written to {args.output_file}")
