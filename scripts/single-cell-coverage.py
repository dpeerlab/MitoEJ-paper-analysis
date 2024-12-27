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
parser.add_argument('--mito_frags_path', type=str, required=True)
parser.add_argument('--barcodes_path', type=str, required=True)
parser.add_argument('--out_dir', type = str, required=True)

args = parser.parse_args()
mito_frags_path = args.mito_frags_path
barcodes_path = args.barcodes_path
out_dir = args.out_dir

len_mt = 16569

mito_genes = pd.read_csv(args.mito_genes_path, index_col=0)
frag = pd.read_csv(f'{mito_frags_path}', index_col=0)

barcodes = pd.read_csv(f'{barcodes_path}', header=None)
barcodes = barcodes[0].tolist()

umis_lst = []
counts_lst= []

c = 0
for barc in barcodes:
    if c % 200 == 0:
        print(c)
        
    b = [barc]
    frag_temp = frag.loc[frag['barcode'].isin(b)]
    cov = pd.DataFrame(list(range(0, len_mt)), columns=['bp'])
    cov['count'] = 0
    cov['umis'] = 0

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

        count = ((frag_temp['start'] <= bp) & (frag_temp['end'] >= bp)).sum()
        umis = frag_temp[((frag_temp['start'] <= bp) & (frag_temp['end'] >= bp))]['umis'].sum()
        cov.loc[idx, 'count'] = count
        cov.loc[idx, 'umis'] = umis
    
    counts_lst.append(cov['count'].values)
    umis_lst.append(cov['umis'].values)
    c += 1

counts_df = pd.DataFrame(counts_lst, index = barcodes)
umis_df = pd.DataFrame(umis_lst, index = barcodes)

counts_df.to_csv(f"{out_dir}counts_df.csv") 
umis_df.to_csv(f"{out_dir}umis_df.csv") 