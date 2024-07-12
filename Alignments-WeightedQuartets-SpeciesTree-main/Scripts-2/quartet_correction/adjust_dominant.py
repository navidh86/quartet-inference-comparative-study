# RUN FROM WITHIN CODE DIRECTORY

from Quartet import Quartet
import random
from copy import deepcopy
import sys

ils = sys.argv[1]
genes = sys.argv[2]
sites = sys.argv[3]
threshold = sys.argv[4] # lower weights get adjusted only
if len(sys.argv) > 5:
    w1 = float(sys.argv[5])
    w2 = float(sys.argv[6])
else:
    w1 = 0.2
    w2 = 0.8


def adjust(est_f, replicate, w1=0.2, w2=0.8): # weighted quartes as arguments
    print("Adjusting replicate", replicate)
    # Read the quartets
    taxa = set()
    covered = {} # quartets that have been covered

    rep = "R" + str(replicate)

    # tq = {}
    # with open(true_f + rep + "/quartets/true.wqrts", "r") as fp:
    #     lines = fp.readlines()
    #     for line in lines:
    #         l = line.split()
    #         temp = Quartet(l[0])
    #         tq[temp] = float(l[1])
    #         # taxa.update(temp.taxa)


    eq = {}
    with open(est_f + rep + "/quartets/est.wqrts", "r") as fp:
        lines = fp.readlines()
        for line in lines:
            l = line.split()
            temp = Quartet(l[0])
            eq[temp] = float(l[1])
            taxa.update(temp.taxa)
            covered[temp] = False

    

    cq = deepcopy(eq)
    dom_q = []
    changes = 0
    for q in cq:
        if covered[q]:
            continue

        q3 = [q]
        q3.extend(q.generate_variants())
        
        best = None
        best_score = -1
        dom = None
        dom_score = -1

        for qq in q3:
            if qq not in eq:
                eq[qq] = 0.0
            covered[qq] = True

            # score = find_score(qq, taxa, eq)
            score = find_score_avg(qq, taxa, eq, w1, w2)
            if score > best_score:
                best_score = score
                best = qq
            
            if eq[qq] > dom_score:
                dom_score = eq[qq]
                dom = qq

        if int(threshold) >= 0 and dom_score > int(threshold) and best != dom:
                changes += 1
                best = dom

        dom_q.append((best.des, eq[best])) # current weight for now

    print("Changes:", changes)

    # write in a file
    output_file = 'est' if int(threshold) < 0 else 'est_'+threshold
    output_file += '_' + str(w1) + '.adqrts'
    with open(est_f + rep + "/quartets/" + output_file, "w") as fp:
        for q in dom_q:
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


def find_score_avg(quartet, taxa, eq, w1=0.2, w2=0.8):
    # w1 = 0.2; w2 = 0.8

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

# if __name__ == "__main__":
#     for r in range(1, 2):
#         print("Adjusting replicate", r)
#         # adjust("15-taxon/100gene-true/", "15-taxon/100gene-100bp/", r)
#         #adjust("37-taxon/noscale.200g.true/", "37-taxon/noscale.100g.500b/", r)
#         # adjust("37-taxon/noscale.200g.true/", "37-taxon/noscale.200g.1500b/", r)
#         # adjust("../avian/avian-1X-truegt/1X-1000-true/", "../avian/avian-1X-estimated-genetrees/1X-1000-1500/", r)
#         adjust("../avian/avian-1X-estimated-genetrees/1X-1000-500/", r)

for r in range(1, 21):
    adjust('../avian/avian-'+ils+'X-estimated-genetrees/'+ils+'X-'+genes+'-'+sites+'/', r, w1, w2)