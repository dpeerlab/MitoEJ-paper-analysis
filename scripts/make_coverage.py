import pandas as pd
import glob
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns
import argparse

custom_params = {'axes.spines.right': False, 'axes.spines.top': False,
                 'figure.dpi':300, 'savefig.dpi':300}
sns.set_theme(style='ticks', rc=custom_params)
sns.set_style({'axes.grid' : False})

parser = argparse.ArgumentParser()
parser.add_argument('--mito_genes_path', type=str, required=True) # mito_genes.csv
parser.add_argument('--mito_fragments_path', type=str, required=True) # mito_fragments.csv
parser.add_argument('--barcodes_path', type=str, required=True) # barcodes.csv
parser.add_argument('--out_path', type=str, required=True)

args = parser.parse_args()

len_mt = 16569
mito_genes = pd.read_csv(args.mito_genes_path, index_col=0)


frag = pd.read_csv(args.mito_fragments_path, index_col=0)
barcodes = pd.read_csv(args.barcodes_path, header=None)
barcodes = barcodes[0].tolist()

print(f'sample: {frag.shape[0]} fragments before filtration')
frag = frag.loc[frag['barcode'].isin(barcodes)]
print(f'sample: {frag.shape[0]} fragments after filtration')

# Create coverage dataframe
cov = pd.DataFrame(list(range(0, len_mt)), columns=['bp'])
cov['count'] = 0

# Determine which gene a base pair is a part of
# chrM if it is not a gene encoding region
cov['gene'] = 'chrM'

for idx, row in mito_genes.iterrows():
    gene = row['gene_name']
    start = row['start']
    end = row['end']
    
    cov.loc[(cov['bp'] >= start) & (cov['bp'] <= end), 'gene'] = gene
    
# Compute the coverage for every base pair
for idx, row in cov.iterrows():
    bp = row['bp']
    
    # See how far into the for loop we are (print base pair)
    if bp % 1000 == 0:
        print(bp)
    
    count = ((frag['start'] <= bp) & (frag['end'] >= bp)).sum()
    cov.loc[idx, 'count'] = count

# Save coverage dataframe to csv
cov.to_csv(args.out_path)