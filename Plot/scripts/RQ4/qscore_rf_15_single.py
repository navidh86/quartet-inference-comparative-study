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
base_folder = '/home/navid/comparative_study/Alignments-WeightedQuartets-SpeciesTree-main'

TAXA = 15

with open(f"{base_folder}/Results/data/Quartet-Scores/{TAXA}-taxon/quartet-scores.pkl", "rb") as fp:
    qscores = pickle.load(fp)

print(qscores)

mapping2 = {
    "astral-regular": "ASTRAL",
    "wQFM-GTF-bucky-boot": "wQFM-GTF-MB",
    "true": "True tree"
}

order = ['wQFM-GTF', 'wQFM-GTF-bucky-boot', 'astral-regular', 'true']
order = [mapping2[x] if x in mapping2 else x for x in order]

REPLICATES = 10 if TAXA == 15 else 20
rf_rates = {}

if TAXA == 15:
    modes = ['100gene-100bp', '1000gene-100bp', '100gene-1000bp', '1000gene-1000bp']
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
    rf_rates[mode]["True tree"] = 0


print(rf_rates)

colors = ['orange', 'blue']
# linestyles = ['-', '--']
linestyles = ['--', '-'] # not used anymore
alphas = [1, 1, 1, 1]
markersize = 60
markers = {
    "wQFM-GTF": 'o',
    "wQFM-GTF-MB": '^',
    "ASTRAL": '*',
    "True tree": 's'
}

capsize = 0.08
label_fontsize = 16
title_fontsize = 18
legend_fontsize = 12
tick_fontsize = 13
ylim_start = 50000
ylim_end = 85000
ytick_lim_start = 50000
ytick_lim_end = 85000 # upto not including
ytick_lim_interval = 2000


modes = ['100gene-100bp', '1000gene-100bp', '100gene-1000bp', '1000gene-1000bp']
### 15 taxon #########################################################
f, ax = plt.subplots(2, 2, figsize=(18, 14), squeeze=True)

for idx, mode in enumerate(modes):
    file_name = f"{TAXA}-{mode}"
    print(f"################ {file_name} ################")

    plt.subplot(2, 2, idx+1)

    x = [rf_rates[mode][method_name] for method_name in order]
    # y_est = [qscores[mode][method_name]['est'] for method_name in order]
    # y_true = [qscores[mode][method_name]['true'] for method_name in order]
    y_est = [qscores[mode][method_name]['est_total'] for method_name in order]
    y_true = [qscores[mode][method_name]['true_total'] for method_name in order]
    x, y_est, y_true = (list(t) for t in zip(*sorted(zip(x, y_est, y_true))))

    print(x)
    print(y_est)
    print(y_true)
    
    plt.plot(x, y_est, color=colors[0], alpha=alphas[idx])
    for method_name in order:
        plt.scatter(rf_rates[mode][method_name], qscores[mode][method_name]['est_total'],
        color=colors[0], marker=markers[method_name], s=markersize)

    
    plt.plot(x, y_true, color=colors[1], alpha=alphas[idx])
    for method_name in order:
        plt.scatter(rf_rates[mode][method_name], qscores[mode][method_name]['true_total'],
        color=colors[1], marker=markers[method_name], s=markersize)

    # plt.xlabel(None)
    # plt.xticks([])

    # plt.yticks(np.arange(ytick_lim_start, ytick_lim_end, ytick_lim_interval), fontsize=tick_fontsize)
    # plt.ylim(ylim_start, ylim_end)

    plt.xlabel("RF Rate", fontsize=label_fontsize)
    plt.ylabel("Quartet Score", fontsize=label_fontsize)
    if idx == 0:
        plt.title("100 Genes", fontsize=title_fontsize)
    elif idx == 1:
        plt.title("1000 Genes", fontsize=title_fontsize)

ax[0][1].yaxis.set_tick_params(labelleft=True)
ax[1][1].yaxis.set_tick_params(labelleft=True)

plt.subplots_adjust(wspace=0.18, hspace=0.12)

# plt.title(f"RF Rate vs Quartet Score ({genes} genes)")

# legend
lines = []
lines.append(Line2D([0,1],[0,1], color=colors[0], label="Estimated gene trees"))    
lines.append(Line2D([0,1],[0,1], color=colors[1], label="True gene trees"))    

# now add markers
lines.append(Line2D([], [], linestyle='', label="")) # for gap
for method_name in order:
    lines.append(Line2D([], [], color='grey', marker=markers[method_name], linestyle='None',
                          markersize=10, label=method_name))

# labels = ['100 bp', '1000 bp', '', 'Estimated gene trees', 'True gene trees', '']
labels = ['Estimated gene trees', 'True gene trees', '']
labels.extend(order)

# bbox to anchor
ba = [0.513, 0.04]

# f.legend(handles=lines, fontsize=legend_fontsize,
#          ncol=7, bbox_to_anchor=ba, loc="lower center")

plt.legend(lines, labels)

out_file = f"{base_folder}/Results/figures/qscore/15taxon-all.pdf"
out_file_png = f"{base_folder}/Results/figures/qscore/15taxon-all.png"
out_file_svg = f"{base_folder}/Results/figures/qscore/15taxon-all.svg"

plt.savefig(out_file, bbox_inches="tight")
plt.savefig(out_file_png, bbox_inches="tight")
plt.savefig(out_file_svg, bbox_inches="tight")

plt.show()
