import os
import shutil
import numpy as np

## GLOBAL ##
TAXON = 15
REPLIACTES = 10
configuration = '100gene-100bp'
input_folder = f'/home/navidh86/Desktop/comparative_study/Alignments-WeightedQuartets-SpeciesTree-main/{TAXON}-taxon/Output/{configuration}'

offset = 300

for mode in ["RAxML"]:
    for replicate_num in range(1, REPLIACTES+1):
        input_file = os.path.join(input_folder, f"R{replicate_num}-{mode}.wqrts")
        output_file = os.path.join(input_folder, f"R{replicate_num}-{mode}.mwqrts")

        with open(input_file, "r") as fp:
            with open(output_file, "w") as wfp:
                lines = fp.readlines()
                for line in lines:
                    line = line.split()
                    qrt, weight = line[0], float(line[1])
                    # print("qrt:", qrt, "weight:", weight)
                    weight += offset
                    wfp.write(qrt + " " + str(weight) + "\n")
                




