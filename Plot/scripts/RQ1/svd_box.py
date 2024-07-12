import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
import xlsxwriter
import sys
from matplotlib.patches import Patch, Rectangle, Polygon
from matplotlib.lines import Line2D
import numpy as np

results = {}
results["method"] = []
results["rf"] = []
results["mode"] = []

mapping2 = {
    'SVD': 'wQFM-SVD (without weights)',
    'wQFM-SVD-exponential': 'wQFM-SVD-exp',
    'wQFM-SVD-reciprocal': 'wQFM-SVD-rec',
    'wQMC-SVD-exponential': 'wQMC-SVD-exp',
    'wQMC-SVD-reciprocal': 'wQMC-SVD-rec'
}

order = ['SVD', 'wQFM-SVD-exponential', 'wQFM-SVD-reciprocal',
         'wQMC-SVD-exponential', 'wQMC-SVD-reciprocal']

order = [mapping2[x] if x in mapping2 else x for x in order]        
hue_order = ['#023858', '#a6bddb', '#74a9cf', '#a1d99b', '#74c476']

if len(sys.argv) > 1:
    TAXA  = int(sys.argv[1])
else:
    TAXA = 15

out_file = "/home/navidh86/Desktop/comparative_study/Alignments-WeightedQuartets-SpeciesTree-main/Results/figures/bar/" + f"RQ1-box-svd-{TAXA}-taxon" + ".pdf"

# SIZES
capsize = 0.08
ylabel_fontsize = 18
title_fontsize = 18
legend_fontsize = 12
tick_fontsize = 13

if TAXA == 10:
    ylim = 0.62
    ytick_lim = 0.7 # upto but not including
elif TAXA == 15:
    ylim = 1
    ytick_lim = 1.1 # upto but not including

if TAXA == 10:
    # 10 lower, 10 higher
    taxon = 10
    for mode in ['lower-ILS', 'higher-ILS']:
        rf_file = f"/home/navidh86/Desktop/comparative_study/Alignments-WeightedQuartets-SpeciesTree-main/Estimated-Species-Trees/RF-rates/{taxon}-taxon/{mode}/RF.txt"

        file_name = f"{taxon}-{mode}"
        print(f"################ {file_name} ################")

        with open(rf_file, "r") as fp:
            for line in fp:
                temp = line[:-1].split()

                method_name = temp[0][temp[0].find('-')+1:] # R1-SVD -> SVD

                if method_name in mapping2:
                    method_name = mapping2[method_name]
                if method_name in order:
                    results["method"].append(method_name)
                    results["rf"].append(float(temp[1]))
                    results["mode"].append(mode)

        data = pd.DataFrame(results)
        # means = data.groupby('method')['rf'].mean().reset_index()

    f, ax = plt.subplots(1, 2, figsize=(20, 8), squeeze=True)

    plt.subplot(1, 2, 1)
    sns.barplot(data=data.query('mode == "lower-ILS" and method in @order'), 
                x='method', y='rf', order=order, width=1, palette=hue_order, capsize=capsize)
    plt.xlabel(None)
    plt.xticks([])
    plt.yticks(np.arange(0.0, ytick_lim, 0.1), fontsize=tick_fontsize)
    plt.ylim(0.0, ylim)
    plt.ylabel("RF Rate", fontsize=ylabel_fontsize)
    plt.title("Lower ILS", fontsize=title_fontsize)


    plt.subplot(1, 2, 2)
    sns.barplot(data=data.query('mode == "higher-ILS" and method in @order'), 
                x='method', y='rf', order=order, width=1, palette=hue_order, capsize=capsize)
    plt.xlabel(None)
    plt.xticks([])
    plt.yticks(np.arange(0.0, ytick_lim, 0.1), fontsize=tick_fontsize)
    plt.ylim(0.0, ylim)
    plt.ylabel("RF Rate", fontsize=ylabel_fontsize)
    plt.title("Higher ILS", fontsize=title_fontsize)

    # show yticks on the right plot too
    ax[1].yaxis.set_tick_params(labelleft=True)

    plt.subplots_adjust(wspace=0.15)

elif TAXA == 15:
    ### 15 taxon #########################################################
    # 100-100....
    taxon = 15
    for mode in ['100gene-100bp', '100gene-1000bp', '1000gene-100bp', '1000gene-1000bp']:
        rf_file = f"/home/navidh86/Desktop/comparative_study/Alignments-WeightedQuartets-SpeciesTree-main/Estimated-Species-Trees/RF-rates/{taxon}-taxon/{mode}/RF.txt"

        file_name = f"{taxon}-{mode}"
        print(f"################ {file_name} ################")
        # out_xl = "/home/navidh86/Desktop/comparative_study/Alignments-WeightedQuartets-SpeciesTree-main/Results/data/RF/" + file_name + ".xlsx    
        # os.makedirs(out_xl.split(f"/{file_name}")[0], exist_ok=True)
        # os.makedirs(out_file.split(f"/{file_name}")[0], exist_ok=True)

        with open(rf_file, "r") as fp:
            for line in fp:
                temp = line[:-1].split()

                method_name = temp[0][temp[0].find('-')+1:] # R1-SVD -> SVD

                if method_name in mapping2:
                    method_name = mapping2[method_name]
                if method_name in order:
                    results["method"].append(method_name)
                    results["rf"].append(float(temp[1]))
                    results["mode"].append(mode)

        data = pd.DataFrame(results)

    # f, ax = plt.subplots(2, 2, sharey='all', figsize=(20, 14), squeeze=True)
    f, ax = plt.subplots(2, 2, figsize=(20, 14), squeeze=True)
    
    plt.subplot(2, 2, 1)
    sns.boxplot(data=data.query('mode == "100gene-100bp" and method in @order'), 
                x='method', y='rf', order=order, width=1, palette=hue_order, orient='v')
    plt.xlabel(None)
    plt.xticks([])
    plt.yticks(np.arange(0.0, ytick_lim, 0.1), fontsize=tick_fontsize)
    plt.ylim(0.0, ylim)
    plt.ylabel("RF Rate", fontsize=ylabel_fontsize)
    plt.title("100 Genes", fontsize=title_fontsize)

    plt.subplot(2, 2, 2)
    sns.boxplot(data=data.query('mode == "1000gene-100bp" and method in @order'), 
                x='method', y='rf', order=order, width=1, palette=hue_order, orient='v')
    plt.xlabel(None)
    plt.xticks([])
    plt.yticks(np.arange(0.0, ytick_lim, 0.1), fontsize=tick_fontsize)
    plt.ylim(0.0, ylim)
    plt.ylabel("RF Rate", fontsize=ylabel_fontsize)
    plt.title("1000 Genes", fontsize=title_fontsize)

    plt.subplot(2, 2, 3)
    sns.boxplot(data=data.query('mode == "100gene-1000bp" and method in @order'), 
                x='method', y='rf', order=order, width=1, palette=hue_order, orient='v')
    plt.xlabel(None)
    plt.xticks([])
    plt.yticks(np.arange(0.0, ytick_lim, 0.1), fontsize=tick_fontsize)
    plt.ylim(0.0, ylim)
    plt.ylabel("RF Rate", fontsize=ylabel_fontsize)
    # plt.title("100 Genes")

    plt.subplot(2, 2, 4)
    sns.boxplot(data=data.query('mode == "1000gene-1000bp" and method in @order'), 
                x='method', y='rf', order=order, width=1, palette=hue_order, orient='v')
    plt.xlabel(None)
    plt.xticks([])
    plt.yticks(np.arange(0.0, ytick_lim, 0.1), fontsize=tick_fontsize)
    plt.ylim(0.0, ylim)
    plt.ylabel("RF Rate", fontsize=ylabel_fontsize)
    # plt.title("100 Genes")

    ax[0][1].yaxis.set_tick_params(labelleft=True)
    ax[1][1].yaxis.set_tick_params(labelleft=True)

    plt.subplots_adjust(wspace=0.15,
                        hspace=0.1)
# plt.subplot_tool()

ncol = 6
if TAXA == 10:
    ba = [0.513, 0.04]
elif TAXA == 15:
    ba = [0.513, 0.07]

legend_elements = []
for idx, item in enumerate(order):
    line = Line2D([0], [0], color=hue_order[idx], marker='s', label=item)
    patch = Rectangle((0, 0), width=0.5, height=2, label=item, color=hue_order[idx],
                    angle=45)
    legend_elements.append(patch)

f.legend(handles=legend_elements, ncol=ncol, loc="lower center", bbox_to_anchor=ba,
            fontsize=legend_fontsize)
plt.savefig(out_file, bbox_inches="tight")
# plt.show()