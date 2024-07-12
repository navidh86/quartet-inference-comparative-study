import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
import xlsxwriter
import sys

results = {}
results["method"] = []
results["rf"] = []

bucky_a = '1-s500'
mapping2 = {
    'wQFM-GTF-dominant-unweighted': 'wQFM-GTF-dom (without weights)',
    'wQMC-GTF-dominant-unweighted': 'wQMC-GTF-dom (without weights)',
    'wQFM-GTF-dominant': 'wQFM-GTF-dom',
    'wQMC-GTF-dominant': 'wQMC-GTF-dom',
    'SVD': 'SVDquartets',
    'wQFM-SVD-exponential': 'wQFM-SVD-exp',
    'wQFM-SVD-reciprocal': 'wQFM-SVD-rec',
    'wQMC-SVD-exponential': 'wQMC-SVD-exp',
    'wQMC-SVD-reciprocal': 'wQMC-SVD-rec',
    f'a{bucky_a}-bucky-mrbayes': 'bucky-MB', 
    f'a{bucky_a}-bucky-RAxML': 'bucky-RAxML',
    f'a{bucky_a}-wQFM-bucky-mrbayes': 'wQFM-bucky-MB',
    f'a{bucky_a}-wQFM-bucky-RAxML': 'wQFM-bucky-RAxML',
    "wQFM-GTF-bucky-boot": "wQFM-GTF-MB",
    "wQMC-GTF-bucky-boot": "wQMC-GTF-MB",
    "wQFM-GTF-boot": "wQFM-GTF-BS",
    "wQMC-GTF-boot": "wQMC-GTF-BS",
}

taxon = 37
# mode = "1X-200-250"
mode = "1X-50-500"
# mode = "1X-200-1000"
# mode = "1000gene-100bp"
# mode = 'higher-ILS'
# mode = "higher-ILS"

if len(sys.argv) > 1:
    taxon = int(sys.argv[1])

if len(sys.argv) > 2:
    mode = sys.argv[2]

print(f"#############{taxon}-{mode}####################")

rf_file = f"/home/navid/comparative_study/Alignments-WeightedQuartets-SpeciesTree-main/Estimated-Species-Trees/RF-rates/{taxon}-taxon/{mode}/RF.txt"
# rf_file = "/home/navidh86/Desktop/comparative_study/Plot/RF-rates/{}-taxon/{}/RF.txt".format(taxon, mode)
# out_file = "/home/navidh86/Desktop/comparative_study/Plot/figures/" + file.split("/")[-1].split(".")[0] + ".png"
file_name = f"{taxon}-{mode}"
# this is default
# out_file = "/home/navidh86/Desktop/comparative_study/Plot/figures/recreated/" + file_name + ".png"
# out_file = "/home/navidh86/Desktop/comparative_study/Plot/figures/recreated/final/" + file_name + ".png" # the last number (if any) is value of w1 * 100
out_xl = "/home/navid/comparative_study/Alignments-WeightedQuartets-SpeciesTree-main/Results/data/RF/" + file_name + ".xlsx"
out_file = "/home/navid/comparative_study/Alignments-WeightedQuartets-SpeciesTree-main/Results/figures/bar/" + file_name + ".png"
os.makedirs(out_xl.split(f"/{file_name}")[0], exist_ok=True)
os.makedirs(out_file.split(f"/{file_name}")[0], exist_ok=True)


with open(rf_file, "r") as fp:
    for line in fp:
        temp = line[:-1].split()

        # results["method"].append(mapping[temp[0].split("/")[0]])
        method_name = temp[0][temp[0].find('-')+1:] # R1-SVD -> SVD
        if method_name in mapping2:
            method_name = mapping2[method_name]

        results["method"].append(method_name)
        results["rf"].append(float(temp[1]))

data = pd.DataFrame(results)
means = data.groupby('method')['rf'].mean().reset_index()
print(means)
# means.reset_index().to_excel(out_xl, engine='xlsxwriter', index=False)
writer = pd.ExcelWriter(out_xl, engine='xlsxwriter')
means.to_excel(writer, 'results', index=False)
means.sort_values(by=['rf']).to_excel(writer, 'sorted', index=False)
data.to_excel(writer, 'All')
writer.save()

exit(0)

# order = ['wQFM-SVD-Rec', 'wQFM-GTF', 'wQMC-SVD-Rec', 'wQMC-GTF']
# order = ['wQFM-SVD-Rec', 'wQFM-SVD-Exp', 'SVD-Tree', 'wQMC-SVD-Exp', 'wQMC-SVD-Rec']
# order = ['wQFM-SVD-Exp', 'wQFM-SVD-Rec', 'wQFM-GTF', 'SVD-Tree', 'wQFM-RAxML', 'wQFM-RAxML-modified', 'wQMC-SVD-Exp', 'wQMC-SVD-Rec', 'wQMC-GTF']
# order = ['wQFM-SVD-Exp', 'wQFM-SVD-Rec', 'wQFM-GTF', 'SVD-Tree', 'wQFM-RAxML', 'wQMC-SVD-Exp', 'wQMC-SVD-Rec', 'wQMC-GTF', 
#          'bucky-pop', 'bucky-conc', 'wQFM-bucky', 'wQMC-bucky',
#          'bucky-pop-boot', 'bucky-conc-boot', 'wQFM-bucky-boot', 'wQMC-bucky-boot']
        #  'bucky-pop-boot-a1', 'wQFM-bucky-boot-a1', 'wQMC-bucky-boot-a1']
# this is default
# order = ['SVD-Tree', 'wQFM-SVD-Exp', 'wQFM-SVD-Rec', 'wQFM-GTF', 'wQFM-RAxML', 'wQMC-SVD-Exp', 'wQMC-SVD-Rec', 'wQMC-GTF']

# 10 taxon
order = ["SVD", "wQFM-SVD-exponential", "wQFM-SVD-reciprocal", "wQFM-GTF", 
         "wQFM-GTF-boot", "wQFM-GTF-adjusted", "wQFM-GTF-boot-adjusted", "wQFM-GTF-distance",
        "wQMC-SVD-exponential", "wQMC-SVD-reciprocal", "wQMC-GTF", "wQMC-GTF-boot", 
        "wQMC-GTF-adjusted", "wQMC-GTF-boot-adjusted", "wQFM-GTF-dominant", 
        "wQFM-GTF-dominant-unweighted", "wQMC-GTF-dominant", "wQMC-GTF-dominant-unweighted",
        "wQMC-GTF-distance"] # "wQFM-RAxML","wQMC-RAxML"
order = [mapping2[x] if x in mapping2 else x for x in order]

bucky_methods = ['bucky-conc-mrbayes', 'bucky-conc-RAxML', 'bucky-mrbayes', 'bucky-RAxML',
                 'wQFM-bucky-mrbayes','wQFM-bucky-RAxML','wQMC-bucky-mrbayes','wQMC-bucky-RAxML'
               ]
# order = []
for method in bucky_methods:
    for alpha in ["1", "1-s500"]:
        order.append(f"a{alpha}-{method}")

# filter out the methods not run yet
order = list(filter(lambda method: method in results["method"], order))

# this is default
# hue_order = ['purple', 'blue', 'green', 'red', 'orange', 'yellow', 'black', 'pink']
# hue_order = ['purple', 'blue', 'green', 'red', 'orange', 'yellow', 'black', 'pink']
hue_order = ['blue', 'green', 'red', 'purple', 'orange', 'black', 'pink', 'cyan', 'yellow', 'grey', 'maroon', 'steelblue']

if mode == "100gene-100bp" and False:
    order.extend(["bucky-mrbayes", "wQFM-bucky-mrbayes", "wQMC-bucky-mrbayes",
				"bucky-RAxML", "wQFM-bucky-RAxML", "wQMC-bucky-RAxML"])
    hue_order.extend(['cyan', 'grey', 'maroon', 'olive', 'brown', 'lime'])

plt.figure(figsize=(12, 8))
plt.gcf().subplots_adjust(bottom=0.15)
# sns.barplot(data=data.query('method in @order'), x='method', y='rf', order=order, palette=hue_order, width=1)
sns.barplot(data=data.query('method in @order'), x='method', y='rf', order=order, width=1)
plt.xticks(rotation=30, size=8)
# plt.legend()
plt.xlabel("Method")
plt.ylabel("RF Rate")
plt.savefig(out_file, bbox_inches="tight")
plt.show()
    