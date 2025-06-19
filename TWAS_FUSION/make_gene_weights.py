import pandas as pd

# --- Load gene lists
with open('gene_list.txt') as f:
    all_genes = set(line.strip() for line in f if line.strip())

with open('unfinished_genes.txt') as f:
    unfinished = set(line.strip() for line in f if line.strip())

genes_to_include = sorted(all_genes - unfinished)

# --- Parse GFF3 for gene coordinates
gff_file = 'Saccharomyces_cerevisiae.R64-1-1.112.gff3'
gene_info = {}
with open(gff_file) as gff:
    for line in gff:
        if line.startswith('#'):
            continue
        cols = line.strip().split('\t')
        if len(cols) < 9: continue
        if cols[2] == 'gene':
            attrs = {kv.split('=')[0]: kv.split('=')[1] for kv in cols[8].split(';') if '=' in kv}
            gene_id = attrs.get('ID', '').replace('gene:', '')
            chrom = cols[0]
            start = int(cols[3])
            end = int(cols[4])
            gene_info[gene_id] = (chrom, start, end)

# --- Write weights_list.txt
out = open('weights_list.txt', 'w')
out.write('WGT\tID\tCHR\tP0\tP1\n')
for gene in genes_to_include:
    if gene not in gene_info:
        print(f'WARNING: {gene} not found in GFF3, skipping.')
        continue
    chrom, start, end = gene_info[gene]
    p0 = max(1, start-200)
    p1 = end+200
    wgt_file = f'{gene}_weights.wgt.RDat'
    out.write(f'{wgt_file}\t{gene}\t{chrom}\t{p0}\t{p1}\n')
out.close()

print('Done! File weights_list.txt created.')
