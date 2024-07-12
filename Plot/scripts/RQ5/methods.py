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

bucky_a = '1-s500'

mapping2 = {
    'SVD': 'SVD-Tree (without weights)',
    "astral": "ASTRAL-old",
    "astral-regular": "ASTRAL",
    f'a{bucky_a}-bucky-mrbayes': 'bucky-mb',
    'wQFM-SVD-reciprocal': 'wQFM-SVD-Rec',
    "wQFM-GTF-bucky-boot": "wQFM-GTF-mb-BS",
    "wQFM-GTF-boot": "wQFM-GTF-RAxML-BS"
}

# order = ['SVD', 'astral', 'astral-regular', f'a{bucky_a}-bucky-mrbayes',
#         'wQFM-SVD-reciprocal', 'wQFM-GTF-bucky-boot', 'wQFM-GTF']

order = ['SVD', 'astral-regular', f'a{bucky_a}-bucky-mrbayes',
        'wQFM-SVD-reciprocal', 'wQFM-GTF-bucky-boot', 'wQFM-GTF',
        'wQFM-GTF-boot']

order = [mapping2[x] if x in mapping2 else x for x in order]        

hue_order = ['#b35806', '#abcdef',
             '#f1a340', '#fee0b6',
             '#d8daeb', '#998ec3',
             '#542788']

# hue_order = ['#b35806', 
#              '#f1a340', '#fee0b6',
#              '#d8daeb', '#998ec3',
#              '#542788']


out_file = "/home/navid/comparative_study/Alignments-WeightedQuartets-SpeciesTree-main/Results/figures/bar/" + f"RQ5-mammalian-methods" + ".pdf"

# SIZES
capsize = 0.08
ylabel_fontsize = 18
title_fontsize = 18
legend_fontsize = 12
tick_fontsize = 13
errorbar="se"

ylim = 0.09
ytick_lim = 0.1

rf_file = f"/home/navid/comparative_study/Alignments-WeightedQuartets-SpeciesTree-main/Estimated-Species-Trees/RF-rates/Mammalian/RF.txt"


with open(rf_file, "r") as fp:
    for line in fp:
        temp = line[:-1].split()

        method_name = temp[0]

        if method_name in mapping2:
            method_name = mapping2[method_name]
        if method_name in order:
            results["method"].append(method_name)
            results["rf"].append(float(temp[1]))

data = pd.DataFrame(results)
# means = data.groupby('method')['rf'].mean().reset_index()

# f, ax = plt.subplots(1, 1, figsize=(26, 8), squeeze=True)

sns.barplot(data=data.query('method in @order'), 
        x='method', y='rf', order=order, width=1, palette=hue_order, capsize=capsize,
        errorbar=errorbar)
plt.xlabel(None)
plt.xticks([])
plt.yticks(np.arange(0.0, ytick_lim, 0.02), fontsize=tick_fontsize)
plt.ylim(0.0, ylim)
plt.ylabel("RF Rate", fontsize=ylabel_fontsize)

plt.title("Mammalian Dataset", fontsize=title_fontsize)


# plt.subplots_adjust(wspace=0.2)
ba = [0.513, 0.04]
ncol = 7

legend_elements = []
for idx, item in enumerate(order):
    # line = Line2D([0], [0], color=hue_order[idx], marker='s', label=item)
    patch = Rectangle((0, 0), width=0.5, height=2, label=item, color=hue_order[idx],
                    angle=45)
    legend_elements.append(patch)

plt.legend(handles=legend_elements, ncol=ncol, loc="lower center", bbox_to_anchor=ba,
            fontsize=legend_fontsize)
plt.savefig(out_file, bbox_inches="tight")
plt.show()