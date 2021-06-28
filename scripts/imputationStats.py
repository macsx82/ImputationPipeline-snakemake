#!/usr/bin/env python
import sys
import matplotlib
import numpy as np
import scipy
import pandas as pd
from matplotlib import pyplot as plt
import argparse
from itertools import cycle
import os


#functions to extract informations on imputation quality
def getMafClass(maf):
    if maf <= 0.01:
        return 'rare'
    elif maf > 0.01 and maf < 0.05:
        return 'low'
    elif maf >= 0.05:
        return 'common'

def getInfoClass(info):
    if info < 0.4:
        return 'lt_04'
    elif info >= 0.04 and info < 0.8:
        return 'gt_04_lt_08'
    elif info >= 0.8:
        return 'gt_08'

################################################################################


#we will read from stdin, but we want to specify file names for table and pdf/figure output
parser = argparse.ArgumentParser()
parser.add_argument('--tab', type=str, help="MANDATORY: Provide full path to the output file for the table report, but only the name prefix.")
parser.add_argument('--fig', type=str, help="MANDATORY: Provide full path to the output file for the plot (PNG) report.")

args=parser.parse_args()
# open files for writing
# out_tab=open(args.tab,'w')
# test='/home/cocca/analyses/imputation/20210613/MOLISANI/06.imputed/MERGED/22/22.tab.INFO'
#we will read from stdin the bcftools formatted file
# chrom_table=pd.read_table(test, sep="\t", names=['CHROM','POS','ID','REF','ALT','AF','INFO_SCORE'])
chrom_table=pd.read_table(sys.stdin, sep="\t", names=['CHROM','POS','ID','REF','ALT','AF','INFO_SCORE'])

# add a MAF column to filter on
chrom_table['MAF'] = [1 - x if x > 0.5 else x for x in chrom_table.AF ]

#now extract some counts
# the INFO score distribution for each chr:

# with all sites
# sites with MAF > 5%
# sites with 1% < MAF <5%
# sites with MAF < 1%

# extract some summary statistics for the INFO score and MAF column
descr_stats = chrom_table[['INFO_SCORE','MAF']].describe()

# get all sites with INFO score eq 0
tot_sites = len(chrom_table)
info_zero = len(chrom_table[chrom_table.INFO_SCORE == 0])

# get all sites with INFO score lt 0.4 and gt 0
info_zero_04 = len(chrom_table[(chrom_table.INFO_SCORE < 0.4) & (chrom_table.INFO_SCORE > 0)])

# get all sites with INFO score gte 0.4
info_04 = len(chrom_table[chrom_table.INFO_SCORE >= 0.4])

# the first table we want to write is this, than
# | CHR | INFO == 0 | 0 < INFO < 0.4 | INFO >=0.4 |
chrom_table['MAF_CLASS'] = chrom_table.MAF.apply(getMafClass)
chrom_table['INFO_CLASS'] = chrom_table.INFO_SCORE.apply(getInfoClass)

# We also want a table with all the numbers
# extract counts for info score stratified by maf classes and info score classes
resume_by_info_maf_classes=chrom_table[['INFO_SCORE','MAF_CLASS','INFO_CLASS']].groupby(['INFO_CLASS','MAF_CLASS']).describe()
# write the table
resume_by_info_maf_classes.to_csv(args.tab + '_by_maf_by_info.csv')
# out_tab.close()

# I want a table with:
# tot site number | tot sites info >=0.4 | ratio of good sites
# overall and by maf class
tot_sites_common = len(chrom_table[(chrom_table.MAF_CLASS == 'common')])
tot_sites_low = len(chrom_table[(chrom_table.MAF_CLASS == 'low')])
tot_sites_rare = len(chrom_table[(chrom_table.MAF_CLASS == 'rare')])
tot_sites_by_classes = [tot_sites, tot_sites_common,tot_sites_low,tot_sites_rare]

info_04_sites_common = len(chrom_table[(chrom_table.MAF_CLASS == 'common') & ((chrom_table.INFO_CLASS == 'gt_04_lt_08') | (chrom_table.INFO_CLASS == 'gt_08')) ])
info_04_sites_low = len(chrom_table[(chrom_table.MAF_CLASS == 'low') & ((chrom_table.INFO_CLASS == 'gt_04_lt_08') | (chrom_table.INFO_CLASS == 'gt_08')) ])
info_04_sites_rare = len(chrom_table[(chrom_table.MAF_CLASS == 'rare') & ((chrom_table.INFO_CLASS == 'gt_04_lt_08') | (chrom_table.INFO_CLASS == 'gt_08')) ])
good_sites_by_classes = [info_04,info_04_sites_common,info_04_sites_low, info_04_sites_rare]


resume = pd.DataFrame(list(zip(tot_sites_by_classes,good_sites_by_classes)), columns=['tot_sites','good_sites'],index=['tot','common','low','rare'])

resume['good_prop'] = (resume.good_sites/resume.tot_sites)*100

resume.to_csv(args.tab + '_by_maf.csv')


# we also want to generate some density plots
# We want to plot the distribution for each MAF slice
# fig, axes = plt.subplots(nrows=3, ncols=1)
# chrom_table.MAF.plot.kde(ax=axes[0])
# chrom_table.INFO_SCORE.plot.kde(ax=axes[1])

# chrom_table[['MAF','INFO_SCORE']].plot.kde(subplots=True)
chrom_table[['MAF','INFO_SCORE']].plot.hist(alpha=0.5, bins=50,subplots=True)
# chrom_table.MAF.plot.hist(alpha=0.5, bins=50,ax=axes[0])
# chrom_table.INFO_SCORE.plot.hist(alpha=0.5, bins=50,ax=axes[1])
plt.legend(loc='best')

# plt.savefig('TEST4.png')
plt.savefig(args.fig)
plt.close()
