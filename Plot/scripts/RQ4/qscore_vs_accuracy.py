import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
import xlsxwriter
import sys
from matplotlib.patches import Patch, Rectangle, Polygon
from matplotlib.lines import Line2D
import numpy as np
import pickle

# base_folder = os.path.dirname(os.path.dirname((os.path.abspath(__file__))))
base_folder = '/home/navidh86/Desktop/comparative_study/Alignments-WeightedQuartets-SpeciesTree-main'

if len(sys.argv) > 1:
    TAXA = int(sys.argv[1])
else:
    TAXA = 15

with open(f"{base_folder}/Results/data/Quartet-Scores/{TAXA}-taxon/quartet-scores.pkl", "rb") as fp:
    qscores = pickle.load(fp)

print(qscores)
exit()

mapping2 = {
    "astral-regular": "ASTRAL",
    "wQFM-GTF-bucky-boot": "wQFM-GTF-mb-BS",
}

order = ['wQFM-GTF', 'wQFM-GTF-bucky-boot', 'astral-regular', 'true']
order = [mapping2[x] if x in mapping2 else x for x in order]

REPLICATES = 10 if TAXA == 15 else 20
rf_rates = {}

if TAXA == 15:
    modes = ['100gene-100bp', '100gene-1000bp', '1000gene-100bp', '1000gene-1000bp']
elif TAXA == 10:
    modes = ['lower-ILS', 'higher-ILS']
elif TAXA == 37:
    # variety = "ILS"
    # if len(sys.argv) > 2:
    #     variety = sys.argv[2]

    # if variety == "ILS":
    #     modes = ['0.5X-200-500', '1X-200-500', '2X-200-500']
    # elif variety == "GENES":
    #     modes = ['1X-100-500', '1X-200-500', '1X-500-500']
    # elif variety == "BP":
    #     modes = ['1X-200-250', '1X-200-500', '1X-200-1000']
    modes = ['0.5X-200-500', '1X-200-500', '2X-200-500', '1X-100-500', '1X-500-500', 
                  '1X-200-250', '1X-200-1000']

for mode in modes:
    rf_rates[mode] = {}
    rf_file = f"{base_folder}/Estimated-Species-Trees/RF-rates/{TAXA}-taxon/{mode}/RF.txt"

    file_name = f"{TAXA}-{mode}"
    print(f"################ {file_name} ################")

    with open(rf_file, "r") as fp:
        for line in fp:
            temp = line[:-1].split()

            method_name = temp[0][temp[0].find('-')+1:] # R1-SVD -> SVD

            if method_name in mapping2:
                method_name = mapping2[method_name]
            if method_name in order:
                if method_name not in rf_rates[mode]:
                    rf_rates[mode][method_name] = float(temp[1])
                else:
                    rf_rates[mode][method_name] += float(temp[1])

    for method_name in rf_rates[mode]:
        rf_rates[mode][method_name] /= REPLICATES
    rf_rates[mode]["true"] = 0


print(rf_rates)


# est plot
if len(sys.argv) > 2:
    mode = sys.argv[2]
else:
    mode = modes[0]

x = [rf_rates[mode][method_name] for method_name in order]
show_total = ""
if len(sys.argv) > 3 and sys.argv[3] == "total":
    y_est = [qscores[mode][method_name]['est_total'] for method_name in order]
    y_true = [qscores[mode][method_name]['true_total'] for method_name in order]
    show_total = "-total"
else:
    y_est = [qscores[mode][method_name]['est'] for method_name in order]
    y_true = [qscores[mode][method_name]['true'] for method_name in order]


# print(x)
# print(y_est)
# print(y_true)
x, y_est, y_true = (list(t) for t in zip(*sorted(zip(x, y_est, y_true))))
# print(x)
# print(y_est)
# print(y_true)

# plt.plot(x, y_est, linestyle='-', marker='o')
plt.plot(y_est, x, linestyle='-', marker='o')
# plt.plot(x, y_true, linestyle='-', marker='o')
plt.plot(y_true, x, linestyle='-', marker='o')

plt.legend(['estimated', 'true'])

# plt.xlabel("RF rate")
# plt.ylabel("Quartet score")
plt.ylabel("RF rate")
plt.xlabel("Quartet score")

plt.title(f"{TAXA} taxa ({mode})")


out_file = f"{base_folder}/Results/figures/qscore/{TAXA}-taxon-{mode}{show_total}.png"

plt.savefig(out_file, bbox_inches="tight")
plt.show()
exit()


results = {}
results["method"] = []
results["rf"] = []
results["mode"] = []


mapping2 = {
    "astral-regular": "ASTRAL",
    "wQFM-GTF-bucky-boot": "wQFM-GTF-mb-BS",
}

order = ['wQFM-GTF', 'wQFM-GTF-bucky-boot', 'astral-regular', 'true']

order = [mapping2[x] if x in mapping2 else x for x in order]        

hue_order = ['#b35806', 
             '#f1a340', '#fee0b6',
             '#d8daeb']

if len(sys.argv) > 1:
    TAXA  = int(sys.argv[1])
else:
    TAXA = 15

if TAXA == 37:
    if len(sys.argv) > 2:
        variety = "-" + sys.argv[2]
    else:
        variety = "-ILS"
else:
    variety = ""


out_file = "/home/navidh86/Desktop/comparative_study/Alignments-WeightedQuartets-SpeciesTree-main/Results/figures/bar/" + f"RQ4-{TAXA}-taxon{variety}" + ".pdf"

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
    ylim = 0.4
    ytick_lim = 0.42 # upto but not including
elif TAXA == 37:
    ylim = 0.09
    ytick_lim = 0.1

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
    sns.barplot(data=data.query('mode == "100gene-100bp" and method in @order'), 
                x='method', y='rf', order=order, width=1, palette=hue_order, capsize=capsize)
    plt.xlabel(None)
    plt.xticks([])
    plt.yticks(np.arange(0.0, ytick_lim, 0.1), fontsize=tick_fontsize)
    plt.ylim(0.0, ylim)
    plt.ylabel("RF Rate", fontsize=ylabel_fontsize)
    plt.title("100 Genes", fontsize=title_fontsize)

    plt.subplot(2, 2, 2)
    sns.barplot(data=data.query('mode == "1000gene-100bp" and method in @order'), 
                x='method', y='rf', order=order, width=1, palette=hue_order, capsize=capsize)
    plt.xlabel(None)
    plt.xticks([])
    plt.yticks(np.arange(0.0, ytick_lim, 0.1), fontsize=tick_fontsize)
    plt.ylim(0.0, ylim)
    plt.ylabel("RF Rate", fontsize=ylabel_fontsize)
    plt.title("1000 Genes", fontsize=title_fontsize)

    plt.subplot(2, 2, 3)
    sns.barplot(data=data.query('mode == "100gene-1000bp" and method in @order'), 
                x='method', y='rf', order=order, width=1, palette=hue_order, capsize=capsize)
    plt.xlabel(None)
    plt.xticks([])
    plt.yticks(np.arange(0.0, ytick_lim, 0.1), fontsize=tick_fontsize)
    plt.ylim(0.0, ylim)
    plt.ylabel("RF Rate", fontsize=ylabel_fontsize)
    # plt.title("100 Genes")

    plt.subplot(2, 2, 4)
    sns.barplot(data=data.query('mode == "1000gene-1000bp" and method in @order'), 
                x='method', y='rf', order=order, width=1, palette=hue_order, capsize=capsize)
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
elif TAXA == 37:
    taxon = 37
    if variety == "-ILS":
        modes = ['0.5X-200-500', '1X-200-500', '2X-200-500'] # ILS
    elif variety == "-BP":
        modes = ['1X-200-250', '1X-200-500', '1X-200-1000'] # BP
    for mode in modes:
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

    f, ax = plt.subplots(1, 3, figsize=(26, 8), squeeze=True)

    for i in range(len(modes)):
        plt.subplot(1, 3, i+1)
        mode = modes[i]
        sns.barplot(data=data.query('mode == @mode and method in @order'), 
                x='method', y='rf', order=order, width=1, palette=hue_order, capsize=capsize)
        plt.xlabel(None)
        plt.xticks([])
        plt.yticks(np.arange(0.0, ytick_lim, 0.02), fontsize=tick_fontsize)
        plt.ylim(0.0, ylim)
        plt.ylabel("RF Rate", fontsize=ylabel_fontsize)
        if variety == "-ILS":
            plt.title(modes[i].split('-')[0], fontsize=title_fontsize)
        elif variety == "-BP":
            plt.title(modes[i].split('-')[2] + " bp", fontsize=title_fontsize)
        if i > 0:
            ax[i].yaxis.set_tick_params(labelleft=True)

    plt.subplots_adjust(wspace=0.2)

ncol = len(order)
if TAXA == 10 or TAXA == 37:
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
plt.show()