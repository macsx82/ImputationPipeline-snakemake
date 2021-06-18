#!/usr/bin/env python
import sys
import matplotlib
import numpy as np
import scipy
import pandas as pd
from matplotlib import pyplot as plt
# import matplotlib.pyplot as plt
import argparse
from itertools import cycle
from natsort import natsorted
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

# #function to create a manhattan plot, code from Stephane Burgeoise
# def plotManhattan(dataframe, **params):
#     # --- Load matplotlib backends ---
# #
#     # get dataframe col names for consistency check
#     header=list(dataframe.columns)
#     regp=dataframe
#     # --- working directory and default parameters ---
#     if 'pvalcolname' in params:
#         pvalcolname = params['pvalcolname']
#     else:
#         pvalcolname = "frequentist_add_pvalue"
# #
#     if 'chromcolname' in params:
#         chromcolname = params['chromcolname']
#     else:
#         chromcolname = "chromosome"
# #
#     if 'poscolname' in params:
#         poscolname = params['poscolname']
#     else:
#         poscolname = "position"
# #    
#     # --- load files ---
#     print("Loading input file(s)...")
#     cols = [chromcolname, poscolname, pvalcolname]
# #   
#     # --- Set plotting defaults ---
#     if 'threshold' in params:
#         threshold = params['threshold']
#     else:
#         threshold = 5E-08
# #
#     if 'colors' in params:
#         c = params['colors'].split(',')
#         if c[0]=="C":
#             from matplotlib.colors import hex2color
#             hexc = ["#000080","#306EFF","#F75D59","#990012"]
#             cs = [hex2color(hexc[i]) for i in range(len(hexc))]
#             colors = cycle(cs)
#         else:
#             colors = cycle(c)
#     else:
#         c = ['black', 'grey']
#         colors = cycle(c)
# #    
#     # --- set ymax ---
#     if 'line' in params:
#         pmin = min(regp[pvalcolname])
#         allmin = min([pmin, threshold])
#         ymax = -np.log10(allmin) + 0.5
#     else:
#         pmin = min(regp[pvalcolname])
#         ymax = -np.log10(pmin) + 0.5
# #
#     if 'ymax' in params:
#         if params['ymax'] < 1:
#             sys.exit("Invalid Ymax value (i.e. entered value < 1)")
# #        
#         if params['ymax'] < ymax:
#             print("Your selected Ymax value is smaller than the actual Ymax and will be ignored")
#         else:
#             ymax = params['ymax']
# #
#     print("Ymax is " + str(ymax))
#     # --- Chromosomes ---
#     #regp[chromcolname] = regp[chromcolname].astype(str)
#     if not regp.loc[regp[chromcolname]=='X',:].empty:
#         regp.loc[regp[chromcolname]=='X',chromcolname] = "23"
# #
#     if not regp.loc[regp[chromcolname]=='Y',:].empty:
#         regp.loc[regp[chromcolname]=='Y',chromcolname] = "24"
# #
#     if not regp.loc[regp[chromcolname]=='XY',:].empty:
#         regp.loc[regp[chromcolname]=='XY',chromcolname] = "25"
# #
#     if not regp.loc[regp[chromcolname]=='MT',:].empty:
#         regp.loc[regp[chromcolname]=='MT',chromcolname] = "26"
# #
#     regp[chromcolname] = regp[chromcolname].astype(int)
# #
#     allc = natsorted(list(set(regp[chromcolname])))
#     chromdisplay = []
#     for x in allc:
#         if x < 23:
#             chromdisplay.append(str(x))
#         else:
#             if x==23:
#                 chromdisplay.append("X")
#             elif x==24:
#                 chromdisplay.append("Y")
#             elif x==25:
#                 chromdisplay.append("XY")
#             elif x==26:
#                 chromdisplay.append("MT")
#             else:
#                 sys.exit("Abnormal chromosome number found: " + str(x))
# #                
#     # --- set x values per chromosome ---
#     xmax = 0
#     cs = {}
#     for chr in allc:
#         colorc = next(colors)
#         xstartc = xmax + 10000000
#         minc = min(regp.loc[regp[chromcolname]==int(chr), poscolname])
#         maxc = max(regp.loc[regp[chromcolname]==int(chr), poscolname])
#         sizec = maxc - minc
#         xmax = xstartc + sizec
#         xtickc = xmax - (sizec/2)
#         cs[str(chr)] = [minc, maxc, sizec, xstartc, xtickc, colorc]
# #
#     xmax = xmax + 10000000
#     # --- Open figure and create axes---
#     # if outformat == ".pdf":
#     #     pp = PdfPages(outfile)
#     #     f = plt.figure()
#     # else:
#     #     f = plt.figure(figsize=dim, dpi=res)
#     f = plt.figure()
#     rect = 0.1,0.1,0.8,0.8
#     ax = f.add_axes(rect)
#     ax.set_ylabel('-log10(p-value)')
#     ax.set_xlabel('chromosome')
#     plt.ylim(0, ymax)
#     plt.xlim(0, xmax)
#     plt.xticks([cs[str(i)][4] for i in allc], chromdisplay, size=4)
#     ax.get_yaxis().tick_left()
#     plt.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
# #
#     if 'line' in params:
#         plt.axhline(y=-math.log10(threshold), color='r', linestyle='--')
# #
#     # --- Add title ---
#     if 'title' in params:
#         plt.title(params['title'])
# #
#     # --- Plot ---
#     for chr in allc:
#         print("Plotting chromosome " + str(chr))
#         xs = np.array(regp.loc[(regp[chromcolname]==int(chr)) & (regp[pvalcolname]>threshold), poscolname]) - cs[str(chr)][0] + cs[str(chr)][3]
#         ys = -np.log10(np.array(regp.loc[(regp[chromcolname]==int(chr)) & (regp[pvalcolname]>threshold), pvalcolname]))
#         ax.plot(np.squeeze(xs), np.squeeze(ys), linestyle='', marker='o', color=cs[str(chr)][5], markerfacecolor=cs[str(chr)][5], markersize=2, markeredgecolor=cs[str(chr)][5], alpha=0.5)
#         # Plot hits below threshold in green
#         sig = regp.loc[(regp[chromcolname]==int(chr)) & (regp[pvalcolname]<=threshold), poscolname]
#         if len(sig) > 0:
#             xs2 = np.array(sig) - cs[str(chr)][0] + cs[str(chr)][3]
#             ys2 = -np.log10(np.array(regp.loc[(regp[chromcolname]==int(chr)) & (regp[pvalcolname]<=threshold), pvalcolname]))
#             ax.plot(np.squeeze(xs2), np.squeeze(ys2), linestyle='', marker='^', color='green', markerfacecolor='green', markersize=3, markeredgecolor='green', alpha=0.5)

    # # --- Save file ---
    # print("Done! Now saving figure...")
    # if outformat == ".pdf":
    #     pp.savefig()
    #     pp.close()
    #     plt.close()
    # elif outformat == ".eps":
    #     plt.savefig(outfile, dpi=res, format='eps', pad_inches=0.1)
    #     plt.close()
    # elif outformat == ".svg":
    #     plt.savefig(outfile, dpi=res, format='svg', pad_inches=0.1)
    #     plt.close()
    # else:
    #     plt.savefig(outfile, dpi=res, format='png', pad_inches=0.1)
    #     plt.close()

# print("Successfully finished. The run took: ",(datetime.now()-startTime))

################################################################################


#we will read from stdin, but we want to specify file names for table and pdf/figure output
parser = argparse.ArgumentParser()
parser.add_argument('--tab', type=str, help="MANDATORY: Provide full path to the output file for the table report.")
parser.add_argument('--fig', type=str, help="MANDATORY: Provide full path to the output file for the plot (PDF) report.")

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
resume_by_info_maf_classes.to_csv(args.tab)
# out_tab.close()
# we also want a manhattan-like plot
# chrom_table['INFO_LOGP'] = (numpy.log10(chrom_table.INFO_SCORE+1))
# chrom_table['INFO_LOGP'] = -(numpy.log10(chrom_table.INFO_SCORE))

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

# plt.savefig('TEST4.pdf')
plt.savefig(args.fig)
plt.close()
# plt.close(fig)


# params={}
# params['chromcolname'] = 'CHROM'
# params['poscolname'] = 'POS'
# params['pvalcolname'] = 'INFO_SCORE'
# # params['mafcolname'] = 'MAF'
# # params['infocolname'] = 'INFO_SCORE'
# params['title'] = 'Info score Manatthan plot'
# params['threshold'] = 0.01
# params['res'] = 300
