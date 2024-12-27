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
parser.add_argument('--sample', type=str, required=True) # mito_genes.csv
parser.add_argument('--mito_genes_path', type=str, required=True)
parser.add_argument('--barcodes_path', type=str, required=True) # barcodes.csv
parser.add_argument('--cov_path', type=str, required=True)
parser.add_argument('--out_path', type=str, required=True)

args = parser.parse_args()
s = args.sample
len_mt = 16569

mito_genes = pd.read_csv(args.mito_genes_path, index_col=0)

genes = ['chrM'] + mito_genes['gene_name'].tolist()
gene_palette = {}
gene_palette['chrM'] = 'gray'
for i, gene in enumerate(genes):
    if gene == 'chrM':
        continue
    else:
        gene_palette[gene] = sns.color_palette("Paired", 14)[i]

barcodes = pd.read_csv(args.barcodes_path, header=None)

cov = pd.read_csv(args.cov_path)
cov['count'] = cov['count']/barcodes.shape[0]

fig, ax = plt.subplots(figsize=(30, 12))

sns.histplot(x=cov['bp'], weights=cov['count'], hue=cov['gene'], bins=len_mt, palette=gene_palette)
sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
ax.set(title=f'Mito coverage from fragments filtered cell barcodes ({s})')
ax.set(xlabel='Position (mitochondrial base pair)', ylabel='Average Coverage per Cell', ylim=(0, None))
plt.savefig(args.out_path, dpi = 300, bbox_inches='tight')