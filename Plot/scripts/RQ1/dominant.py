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
    'wQFM-GTF-dominant-unweighted': 'wQFM-GTF-dom (without weights)',
    'wQMC-GTF-dominant-unweighted': 'wQMC-GTF-dom (without weights)',
    'wQFM-GTF-dominant': 'wQFM-GTF-dom',
    'wQMC-GTF-dominant': 'wQMC-GTF-dom',
    'wQFM-GTF': 'wQFM-GTF-all',
    'wQMC-GTF': 'wQMC-GTF-all', # all will be used later on, renamed to wQFM-GTF
}

order = ['wQFM-GTF-dominant', 'wQFM-GTF-dominant-unweighted', 'wQFM-GTF',
        'wQMC-GTF-dominant', 'wQMC-GTF-dominant-unweighted', 'wQMC-GTF']
order = [mapping2[x] if x in mapping2 else x for x in order]        
hue_order = ['#deebf7', '#9ecae1', '#3182bd', '#e5f5e0', '#a1d99b', '#31a354']

if len(sys.argv) > 1:
    TAXA  = int(sys.argv[1])
else:
    TAXA = 15

if TAXA == 37:
    # BP-new for including 50 base pair config as well
    if len(sys.argv) > 2:
        variety = "-" + sys.argv[2]
    else:
        variety = "-ILS"
else:
    variety = ""

base_folder = "/home/navid/comparative_study/Alignments-WeightedQuartets-SpeciesTree-main"
out_dir = f"{base_folder}/Results/figures/bar"
os.makedirs(out_dir, exist_ok=True)
out_file = f"{out_dir}/RQ1a-dom-{TAXA}-taxon{variety}" + ".pdf"
out_file_png = f"{out_dir}/RQ1a-dom-{TAXA}-taxon{variety}" + ".png"
out_file_svg = f"{out_dir}/RQ1a-dom-{TAXA}-taxon{variety}" + ".svg"

# SIZES
capsize = 0.08
ylabel_fontsize = 18
title_fontsize = 18
legend_fontsize = 12
tick_fontsize = 13
errorbar = "se" # ("ci", 95)

if TAXA == 10:
    ylim = 0.62
    ytick_lim = 0.7 # upto but not including
elif TAXA == 15:
    ylim = 0.47
    ytick_lim = 0.6 # upto but not including
elif TAXA == 37:
    if variety == "-GENES":
        ylim = 0.102
        ytick_lim = 0.11
    elif variety == "-BP-new":
        ylim = 0.245
        ytick_lim = 0.25
    elif variety == "-ILS":
        ylim = 0.102
        ytick_lim = 0.12
    else:
        ylim = 0.11
        ytick_lim = 0.12

if TAXA == 10:
    # 10 lower, 10 higher
    taxon = 10
    for mode in ['lower-ILS', 'higher-ILS']:
        rf_file = f"{base_folder}/Estimated-Species-Trees/RF-rates/{taxon}-taxon/{mode}/RF.txt"

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
                x='method', y='rf', order=order, width=1, palette=hue_order, capsize=capsize, errorbar=errorbar)
    plt.xlabel(None)
    plt.xticks([])
    plt.yticks(np.arange(0.0, ytick_lim, 0.1), fontsize=tick_fontsize)
    plt.ylim(0.0, ylim)
    plt.ylabel("RF Rate", fontsize=ylabel_fontsize)
    plt.title("Lower ILS", fontsize=title_fontsize)


    plt.subplot(1, 2, 2)
    sns.barplot(data=data.query('mode == "higher-ILS" and method in @order'), 
                x='method', y='rf', order=order, width=1, palette=hue_order, capsize=capsize, errorbar=errorbar)
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
    taxon = 15
    for mode in ['100gene-100bp', '100gene-1000bp', '1000gene-100bp', '1000gene-1000bp']:
        rf_file = f"{base_folder}/Estimated-Species-Trees/RF-rates/{taxon}-taxon/{mode}/RF.txt"

        file_name = f"{taxon}-{mode}"
        print(f"################ {file_name} ################")
        # out_xl = "/home/navid/comparative_study/Alignments-WeightedQuartets-SpeciesTree-main/Results/data/RF/" + file_name + ".xlsx    
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
    sns.barplot(data=data.query('mode == "100gene-100bp" and method in @order'), 
                x='method', y='rf', order=order, width=1, palette=hue_order, capsize=capsize, errorbar=errorbar)
    plt.xlabel(None)
    plt.xticks([])
    plt.yticks(np.arange(0.0, ytick_lim, 0.1), fontsize=tick_fontsize)
    plt.ylim(0.0, ylim)
    plt.ylabel("RF Rate", fontsize=ylabel_fontsize)
    plt.title("100 Genes", fontsize=title_fontsize)

    plt.subplot(2, 2, 2)
    sns.barplot(data=data.query('mode == "1000gene-100bp" and method in @order'), 
                x='method', y='rf', order=order, width=1, palette=hue_order, capsize=capsize, errorbar=errorbar)
    plt.xlabel(None)
    plt.xticks([])
    plt.yticks(np.arange(0.0, ytick_lim, 0.1), fontsize=tick_fontsize)
    plt.ylim(0.0, ylim)
    plt.ylabel("RF Rate", fontsize=ylabel_fontsize)
    plt.title("1000 Genes", fontsize=title_fontsize)

    plt.subplot(2, 2, 3)
    sns.barplot(data=data.query('mode == "100gene-1000bp" and method in @order'), 
                x='method', y='rf', order=order, width=1, palette=hue_order, capsize=capsize, errorbar=errorbar)
    plt.xlabel(None)
    plt.xticks([])
    plt.yticks(np.arange(0.0, ytick_lim, 0.1), fontsize=tick_fontsize)
    plt.ylim(0.0, ylim)
    plt.ylabel("RF Rate", fontsize=ylabel_fontsize)
    # plt.title("100 Genes")

    plt.subplot(2, 2, 4)
    sns.barplot(data=data.query('mode == "1000gene-1000bp" and method in @order'), 
                x='method', y='rf', order=order, width=1, palette=hue_order, capsize=capsize, errorbar=errorbar)
    plt.xlabel(None)
    plt.xticks([])
    plt.yticks(np.arange(0.0, ytick_lim, 0.1), fontsize=tick_fontsize)
    plt.ylim(0.0, ylim)
    plt.ylabel("RF Rate", fontsize=ylabel_fontsize)
    # plt.title("100 Genes")

    ax[0][1].yaxis.set_tick_params(labelleft=True)
    ax[1][1].yaxis.set_tick_params(labelleft=True)

    plt.subplots_adjust(wspace=0.15, hspace=0.1)

elif TAXA == 37:
    taxon = 37
    if variety == "-ILS":
        modes = ['0.5X-200-500', '1X-200-500', '2X-200-500'] # ILS
    elif variety == "-BP":
        modes = ['1X-200-250', '1X-200-500', '1X-200-1000'] # BP
    elif variety == "-BP-new":
        modes = ['1X-200-50', '1X-200-250', '1X-200-500', '1X-200-1000'] # BP with 50 bp included
    elif variety == "-GENES":
        modes = ['1X-100-500', '1X-200-500', '1X-500-500'] # GENES
        # modes = ['1X-50-500', '1X-100-500', '1X-200-500'] # GENES
    for mode in modes:
        rf_file = f"{base_folder}/Estimated-Species-Trees/RF-rates/{taxon}-taxon/{mode}/RF.txt"

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

    f, ax = plt.subplots(1, len(modes), figsize=(26, 8), squeeze=True)

    for i in range(len(modes)):
        plt.subplot(1, len(modes), i+1)
        mode = modes[i]
        sns.barplot(data=data.query('mode == @mode and method in @order'), 
                x='method', y='rf', order=order, width=1, palette=hue_order, capsize=capsize, errorbar=errorbar)
        plt.xlabel(None)
        plt.xticks([])
        plt.yticks(np.arange(0.0, ytick_lim, 0.02), fontsize=tick_fontsize)
        plt.ylim(0.0, ylim)
        plt.ylabel("RF Rate", fontsize=ylabel_fontsize)
        if variety == "-ILS":
            plt.title(modes[i].split('-')[0], fontsize=title_fontsize)
        elif variety == "-BP" or variety == "-BP-new":
            plt.title(modes[i].split('-')[2] + " bp", fontsize=title_fontsize)
        elif variety == "-GENES":
            plt.title(modes[i].split('-')[1] + " genes", fontsize=title_fontsize)
        if i > 0:
            ax[i].yaxis.set_tick_params(labelleft=True)

    plt.subplots_adjust(wspace=0.2)

ncol = 6
if TAXA == 10 or TAXA == 37:
    ba = [0.513, 0.04]
elif TAXA == 15:
    ba = [0.513, 0.07]

legend_elements = []
for idx, item in enumerate(order):
    # line = Line2D([0], [0], color=hue_order[idx], marker='s', label=item)
    patch = Rectangle((0, 0), width=0.5, height=2, label=item, color=hue_order[idx], angle=45)
    legend_elements.append(patch)

f.legend(handles=legend_elements, ncol=ncol, loc="lower center", bbox_to_anchor=ba, fontsize=legend_fontsize)
plt.savefig(out_file, bbox_inches="tight")
plt.savefig(out_file_png, bbox_inches="tight")
plt.savefig(out_file_svg, bbox_inches="tight")
# plt.show()