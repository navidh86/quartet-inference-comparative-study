# RUN FROM WITHIN CODE DIRECTORY

from Quartet import Quartet
import random
from copy import deepcopy
import sys

# ils = sys.argv[1]
# genes = sys.argv[2]
# sites = sys.argv[3]
# threshold = sys.argv[4] # lower weights get adjusted only
# if len(sys.argv) > 5:
#     w1 = float(sys.argv[5])
#     w2 = float(sys.argv[6])
# else:
#     w1 = 0.8
#     w2 = 0.2

input_file = sys.argv[1]
output_file = sys.argv[2]
if len(sys.argv) > 3:
    w1 = float(sys.argv[3])
    w2 = float(sys.argv[4])
else:
    w1 = 0.8
    w2 = 0.2


def adjust():
    print(f"Adjusting {input_file}")

    taxa = set()
    covered = {} # quartets that have been covered

    eq = {} # input quartets
    with open(input_file, "r") as fp:
        lines = fp.readlines()
        for line in lines:
            l = line.split()
            temp = Quartet(l[0])
            eq[temp] = float(l[1])
            taxa.update(temp.taxa)
            covered[temp] = False

    cq = deepcopy(eq)
    dom_q = []
    weights = [] # list of all quartets with their weights
    wq = {} # used to hold newly computed weights
    changes = 0

    for q in cq:
        if covered[q]:
            continue

        q3 = [q]
        q3.extend(q.generate_variants())

        total_score = 0
        total_quartets = 0
        for qq in q3:
            if qq not in eq:
                eq[qq] = 0.0
            covered[qq] = True
            total_quartets += eq[qq]

            # score = find_score(qq, taxa, eq)
            wq[qq] = find_score_avg(qq, taxa, eq, w1, w2)
            total_score += wq[qq]

        for qq in q3:
            weights.append((qq.des, wq[qq]/total_score*total_quartets))

    # print("Changes:", changes)
    weights.sort()
    # write in a file
    with open(output_file, "w") as fp:
        for q in weights:
            fp.write(q[0] + " " + str(q[1]) + "\n")


def find_score(quartet, taxa, eq):
    score = 0
    new_taxa = set()
    for taxon in taxa:
        if taxon not in quartet.taxa:
            new_taxa.add(taxon)

    for taxon in new_taxa:

        temp = Quartet("((" + taxon + "," + quartet.taxa[1] + "),(" + quartet.taxa[2] + "," + quartet.taxa[3] + "));")
        if temp in eq:
            score += eq[temp]
        temp = Quartet("((" + quartet.taxa[0] + "," + taxon + "),(" + quartet.taxa[2] + "," + quartet.taxa[3] + "));")
        if temp in eq:
            score += eq[temp]
        temp = Quartet("((" + quartet.taxa[0] + "," + quartet.taxa[1] + "),(" + taxon + "," + quartet.taxa[3] + "));")
        if temp in eq:
            score += eq[temp]
        temp = Quartet("((" + quartet.taxa[0] + "," + quartet.taxa[1] + "),(" + quartet.taxa[2] + "," + taxon + "));")
        if temp in eq:
            score += eq[temp]        

    return score


def find_score_avg(quartet, taxa, eq, w1=0, w2=1):
    score = 0
    new_taxa = set()
    for taxon in taxa:
        if taxon not in quartet.taxa:
            new_taxa.add(taxon)

    for taxon in new_taxa:
        temp = Quartet("((" + taxon + "," + quartet.taxa[1] + "),(" + quartet.taxa[2] + "," + quartet.taxa[3] + "));")
        if temp in eq:
            score += eq[temp]
        temp = Quartet("((" + quartet.taxa[0] + "," + taxon + "),(" + quartet.taxa[2] + "," + quartet.taxa[3] + "));")
        if temp in eq:
            score += eq[temp]
        temp = Quartet("((" + quartet.taxa[0] + "," + quartet.taxa[1] + "),(" + taxon + "," + quartet.taxa[3] + "));")
        if temp in eq:
            score += eq[temp]
        temp = Quartet("((" + quartet.taxa[0] + "," + quartet.taxa[1] + "),(" + quartet.taxa[2] + "," + taxon + "));")
        if temp in eq:
            score += eq[temp]        



    return w1 * eq[quartet] + w2 * score / (4 * len(new_taxa))

adjust()