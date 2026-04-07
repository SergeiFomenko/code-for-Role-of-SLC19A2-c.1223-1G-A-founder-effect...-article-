from cyvcf2 import VCF
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as patches
from matplotlib.lines import Line2D
parser = argparse.ArgumentParser(
    description="""Haplotype analysis""",
    formatter_class=argparse.RawTextHelpFormatter
)
parser.add_argument(
    '--multi_vcf',
    type=Path,
    required=True,
    help="Common vcf with heterogyous carriers and homozygous samples"
)
parser.add_argument(
    '--nokin_subset',
    type=Path,
    required=True,
    help="Subset of non-related individuals from cohort"
)
parser.add_argument(
    '--variant',
    type=str,
    required=True,
    help="Variant data in format CHROM_POS_REF_ALT"
)
parser.add_argument(
    '--region',
    type=str,
    required=True,
    help="Region of homozygosity in format CHROM_START_END"
)
parser.add_argument(
    '--results_path',
    type=Path,
    required=True,
    help="Resulting stats path"
)
parser.add_argument(
    '--plot_path',
    type=Path,
    required=True,
    help="Resulting plot path"
)
parser.add_argument(
    '--genes',
    type=Path,
    required=True,
    help="Path to genes for plot annotation"
)
args = parser.parse_args()
multi_vcf = VCF(args.multi_vcf)
pop_samples = multi_vcf.samples
region = args.region.split('_')
if len(region) != 3:
    raise ValueError("Region must be CHROM_START_END")
region[1] = int(region[1])
region[2] = int(region[2])
padded_region = (max(region[1] - 500000,0),region[2] + 500000)
variant = args.variant.split('_')
if len(variant) != 4:
    raise ValueError("Variant must be CHROM_POS_REF_ALT")
variant[1] = int(variant[1])
genes = args.genes
results_path = args.results_path
plot_path = args.plot_path
all_var = []
for v in multi_vcf(f'{region[0]}:{padded_region[0]}-{padded_region[1]}'):
    all_var.append(['_'.join([v.CHROM,str(v.POS),v.REF,v.ALT[0]])] + list(v.gt_types))
multi_vcf.close()
all_var = pd.DataFrame(all_var, columns=['var_id'] + pop_samples)
ind29_dict = {}
with open(args.nokin_subset,'r') as ind29:
    ind29_dict = {line.strip():True for line in ind29.readlines()}
var_gts = all_var[all_var['var_id'] == '_'.join([str(val) for val in variant])].reset_index(drop=True).T
hetero_carriers = var_gts[var_gts[0] == 1].index.to_list()
homoref_carriers = var_gts[var_gts[0] == 0].index.to_list()
homo_patients = var_gts[var_gts[0] == 3].index.to_list()
haplotype_homo_vars = all_var[(all_var[homo_patients] == 3).all(axis=1) | (all_var[homo_patients] == 0).all(axis=1)]
haplotype_homo_vars['POS'] = haplotype_homo_vars['var_id'].str.split('_',expand=True)[1].astype(int)
snps = all_var[(all_var['var_id'].str.split('_',expand=True)[2].str.len() == 1) &\
    (all_var['var_id'].str.split('_',expand=True)[3].str.len() == 1)]
snps['patient_heterorate'] = ((snps[homo_patients] == 1) * 1).mean(axis=1)
snps['POS'] = snps['var_id'].str.split('_',expand=True)[1].astype(int)
alt_homovars = snps[(snps[homo_patients] == 3).all(axis=1)]
alt_homovars['mismatch_carriers'] = (alt_homovars[[zlims for zlims in hetero_carriers if ind29_dict.get(zlims.lstrip('0'))]] == 0).mean(axis=1)
alt_homovars['mismatch_non-carriers'] = (alt_homovars[[zlims for zlims in homoref_carriers if ind29_dict.get(zlims.lstrip('0'))]] == 0).mean(axis=1)
alt_homovars['haplotype_gt'] = '1/1'
ref_homovars = snps[(snps[homo_patients] == 0).all(axis=1)]
ref_homovars['mismatch_carriers'] = (ref_homovars[[zlims for zlims in hetero_carriers if ind29_dict.get(zlims.lstrip('0'))]] == 3).mean(axis=1)
ref_homovars['mismatch_non-carriers'] = (ref_homovars[[zlims for zlims in homoref_carriers if ind29_dict.get(zlims.lstrip('0'))]] == 3).mean(axis=1)
ref_homovars['haplotype_gt'] = '0/0'
mismatch_rate = pd.concat([alt_homovars[['var_id','mismatch_carriers','mismatch_non-carriers','haplotype_gt']],
           ref_homovars[['var_id','mismatch_carriers','mismatch_non-carriers','haplotype_gt']]
          ])
mismatch_rate['POS'] = mismatch_rate['var_id'].str.split('_',expand=True)[1].astype(int)
mismatch_rate = mismatch_rate.sort_values(by='POS')
mismatch_rate['var_gt'] = mismatch_rate['var_id'] + ' ' + mismatch_rate['haplotype_gt']
snps = snps[((snps[snps.columns[snps.columns.str.contains('^0')]] == 2).mean(axis=1) < 0.05) &\
    ((snps[homo_patients] == 2).mean(axis=1) < 0.05)]
snps['homoref_rate_patients'] = (snps[homo_patients] == 0).sum(axis=1) / len(homo_patients)
snps['homoalt_rate_patients'] = (snps[homo_patients] == 3).sum(axis=1) / len(homo_patients)
snps['patients_match_percent'] = snps[['homoref_rate_patients','homoalt_rate_patients']].max(axis=1)



genes = pd.read_csv(genes,sep='\t')[['#chrom','chromStart','chromEnd','geneName2']].fillna('')
orange_dot = Line2D(
    [0], [0],
    marker='o',
    color='w',
    markerfacecolor='orange',
    markeredgecolor='black',
    markersize=6,
    linestyle='None'
)
custom_handle = Line2D([0], [0], color='black', linestyle='--', label='Patients haplotype\nchr1:168,981,142–171,250,829')
fig, axs = plt.subplots(4, 1, sharex=True,figsize=(18, 6))
snps = snps[((snps['POS'] >= padded_region[0]) & (snps['POS'] <= padded_region[1]))].sort_values(by='POS')
mismatch_rate = mismatch_rate[((mismatch_rate['POS'] >= region[1]) & (mismatch_rate['POS'] <= region[2]))].sort_values(by='POS')
axs[0].plot(snps['POS'],1- snps['patients_match_percent'], color=sns.husl_palette(n_colors=3)[0])
axs[0].scatter(variant[1], 0.01, zorder=3, color='orange',edgecolor='black')
axs[0].axvline(x=region[1], color='black', linestyle='--')
axs[0].axvline(x=region[2], color='black', linestyle='--')
axs[1].scatter(variant[1], 0.003, zorder=3, color='orange',edgecolor='black')
axs[1].plot(mismatch_rate['POS'],mismatch_rate['mismatch_carriers'], color=sns.husl_palette(n_colors=3)[1])
axs[1].axvline(x=region[1], color='black', linestyle='--')
axs[1].axvline(x=region[2], color='black', linestyle='--')
axs[2].plot(mismatch_rate['POS'],mismatch_rate['mismatch_non-carriers'], color=sns.husl_palette(n_colors=3)[2])
axs[2].scatter(variant[1], 0.05, color='orange',zorder=3,edgecolor='black')
axs[2].axvline(x=region[1], color='black', linestyle='--')
axs[2].axvline(x=region[2], color='black', linestyle='--')
plt.xlabel('GRCh38 chr1 position, megabases',fontsize=14)
step = 500000
plt.xticks([1000000 * round(pos / 1000000,1) for pos in range(snps['POS'].min(), snps['POS'].max(),step)],
           labels=[round(pos / 1000000,1) for pos in range(snps['POS'].min(), snps['POS'].max(),step)])
axs[0].set_ylabel('percent of\nnon-matching genotypes\nbetween c.1223+1G>A patients',rotation=0,fontsize=12)
axs[0].yaxis.set_label_coords(-0.1, 0.25) 
axs[1].set_ylabel('heterozygous c.1223+1G>A\ncarriers mismatch rate\nwith patients haplotype',rotation=0,fontsize=12)
axs[1].yaxis.set_label_coords(-0.1, 0.25)
axs[2].set_ylabel('non-carriers mismatch rate\nwith patients haplotype',rotation=0,fontsize=12)
axs[2].yaxis.set_label_coords(-0.1, 0.25)
axs[3].set_ylabel('Genes in region',rotation=0,fontsize=12)
axs[3].yaxis.set_label_coords(-0.1, 0.25)



# Subset genes to plotted region
genes_sub = genes[
    (genes['chromEnd'] >= padded_region[0]) &
    (genes['chromStart'] <= padded_region[1])
].copy()

# Sort genes by start
genes_sub = genes_sub.sort_values('chromStart')
min_len = 10000
genes_sub['length'] = genes_sub['chromEnd'] - genes_sub['chromStart']
genes_to_label = genes_sub[genes_sub['length'] > min_len]
ax = axs[3]

height = 0.6
tracks = []  # stores last end position for each track
texts = []   # store text objects for adjustText

for _, row in genes_sub.iterrows():
    start = row['chromStart']
    end = row['chromEnd']
    gene_name = row['geneName2']
    
    # --- interval packing (multi-track) ---
    for i, last_end in enumerate(tracks):
        if start > last_end:
            track_id = i
            tracks[i] = end
            break
    else:
        track_id = len(tracks)
        tracks.append(end)
    
    y = track_id

    # --- draw rectangle ---
    rect = patches.Rectangle(
        (start, y),
        end - start,
        height,
        facecolor='lightgrey',
        edgecolor='black'
    )
    ax.add_patch(rect)

    # --- add text (store for adjustment) ---
    txt = ax.text(
        (start + end) / 2,
        y + height + 0.1,
        gene_name,
        ha='center',
        va='bottom',
        fontsize=6,
        rotation=15
    )
    texts.append(txt)



# Formatting
ax.set_ylim(-0.5, len(tracks) + 1.5)
ax.set_yticks([])
# Formatting
#ax.set_ylim(-0.5, y_level + 2)
fig.legend(handles=[orange_dot, custom_handle],labels=['SLC19A2 c.1223+1G>A\nvariant position',
                                                      f'Patients haplotype\n{region[0]}:{region[1]}–{region[2]}'],
           bbox_to_anchor=(1.2, 0.93),fontsize=12)
plt.xlim(padded_region[0], padded_region[1])
plt.suptitle(f"Region {region[0]}:{region[1]}–{region[2]} in SLC19A2 c.1223+1G>A homozygous patients and Ingush cohort", y=1,fontsize=18)
plt.tight_layout()

plt.savefig(plot_path,bbox_inches='tight',dpi=600)
with pd.ExcelWriter(results_path) as writer:
    snps[['var_id','patient_heterorate', 'homoref_rate_patients', 'homoalt_rate_patients', 'patients_match_percent']]\
        .to_excel(writer, sheet_name='snp_concordance_patients', index=False)
    mismatch_rate[['var_id', 'mismatch_carriers', 'mismatch_non-carriers', 'haplotype_gt']]\
        .to_excel(writer, sheet_name='mismatch_rate_population', index=False)

genotype_concordance = snps.loc[((snps['POS'] >= region[1]) & (snps['POS'] <= region[2])),'patients_match_percent'].mean()
nonrel_carriers = [zlims for zlims in hetero_carriers if ind29_dict.get(zlims.lstrip('0'))]
nonrel_homoref = [zlims for zlims in homoref_carriers if ind29_dict.get(zlims.lstrip('0'))]
mismatch_noncarrier = mismatch_rate.loc[((mismatch_rate['POS'] >= region[1]) & (mismatch_rate['POS'] <= region[2])),'mismatch_non-carriers'].mean()
mismatch_carriers = mismatch_rate.loc[((mismatch_rate['POS'] >= region[1]) & (mismatch_rate['POS'] <= region[2])),'mismatch_carriers'].mean()
print(f'Mean genotype concordance for {len(homo_patients)} homozygous carriers of variant {"_".join([str(val) for val in variant])} in region {region[0]}:{region[1]}–{region[2]} is: {genotype_concordance}')
print(f'Mean genotype mismatch rate for {len(nonrel_carriers)} non-related heterozyous carriers of variant {"_".join([str(val) for val in variant])} in region {region[0]}:{region[1]}–{region[2]} is: {mismatch_carriers}')
print(f'Mean genotype mismatch rate for {len(nonrel_homoref)} non-related non-carriers of variant {"_".join([str(val) for val in variant])} in region {region[0]}:{region[1]}–{region[2]} is: {mismatch_noncarriers}')
