import sys
import subprocess as sp
import itertools
import time

start_time = time.time()

print("computing distance based wqrts...")

# DECLARATIONS
if len(sys.argv) < 3:
    exit("Specify input and output file")

input_gt_file = sys.argv[1]
output_quartets_file = sys.argv[2]
# all_gt_file = "true_all_gt.tre"

def get_qstring(quartet):
    return f"(({quartet[0]},{quartet[1]}),({quartet[2]},{quartet[3]}));"

# create the list of taxa first
with open(input_gt_file, "r") as fp:
    trees = fp.readlines()
    if trees[-1] == "":
        trees = trees[:-1]
    trees = [tree[:-1] for tree in trees]

gt_count = len(trees)

sp.check_output("touch tree.temp", shell=True)
cmd = f"echo '{trees[0]}' > tree.temp"
sp.check_output(cmd, shell=True)

cmd = f"nw_distance -n tree.temp | cut -f 1 | sort | paste -sd ' '"
taxa = sp.check_output(cmd, shell=True).decode("utf-8")[:-1].split(" ")

# to store all weights
quartets = {}
distances = []
cnt = 0

# first compute all pairwise distances
for idx, tree in enumerate(trees):
    distances.append({}) # pairwise distance for this tree
    cmd = f"echo '{tree}' > tree.temp"
    sp.check_output(cmd, shell=True)

    for pair in itertools.combinations(taxa, 2):
        cmd = f"nw_distance -m l tree.temp {pair[0]} {pair[1]} | grep -v e | paste -sd+ | bc"
        # dist = float(sp.check_output(cmd, shell=True).decode('utf-8')[:-1])
        try:
            d = sp.check_output(cmd, shell=True).decode('utf-8')[:-1]
            if d == '':
                dist = 0
            else:
                dist = float(d)
        except:
            print(idx)
            print(pair)
            print(d)
            print(tree)
            exit(f"error occurred on tree {idx} for pair {pair}")
        distances[idx][pair] = dist

    # now loop over all quartets
    for quartet in itertools.combinations(taxa, 4):
        # generate three variants
        qq = [0] * 3
        scores = [0] * 3
        qq[0] = [quartet[0], quartet[1], quartet[2], quartet[3]]
        qq[1] = [quartet[0], quartet[2], quartet[1], quartet[3]]
        qq[2] = [quartet[0], quartet[3], quartet[1], quartet[2]]

        q = qq[0]
        scores[0] = distances[idx][(q[0], q[2])] + \
                    distances[idx][(q[0], q[3])] + \
                    distances[idx][(q[1], q[2])] + \
                    distances[idx][(q[1], q[3])] - \
                    distances[idx][(q[0], q[1])] - \
                    distances[idx][(q[2], q[3])]
        
        scores[1] = distances[idx][(q[0], q[1])] + \
                    distances[idx][(q[0], q[3])] + \
                    distances[idx][(q[1], q[2])] + \
                    distances[idx][(q[2], q[3])] - \
                    distances[idx][(q[0], q[2])] - \
                    distances[idx][(q[1], q[3])]
        
        scores[2] = distances[idx][(q[0], q[1])] + \
                    distances[idx][(q[0], q[2])] + \
                    distances[idx][(q[1], q[3])] + \
                    distances[idx][(q[2], q[3])] - \
                    distances[idx][(q[0], q[3])] - \
                    distances[idx][(q[1], q[2])]

        for idx2, q in enumerate(qq):  
            qstring = get_qstring(q)
            if qstring not in quartets:
                quartets[qstring] = 0
            quartets[qstring] += scores[idx2]
            cnt += 1
       
# normalize the weights
for quartet in itertools.combinations(taxa, 4):
    # generate three variants
    qq = [0] * 3
    qq[0] = get_qstring([quartet[0], quartet[1], quartet[2], quartet[3]])
    qq[1] = get_qstring([quartet[0], quartet[2], quartet[1], quartet[3]])
    qq[2] = get_qstring([quartet[0], quartet[3], quartet[1], quartet[2]])

    total = sum([quartets[q] for q in qq])
    for q in qq:
        quartets[q] = quartets[q] * gt_count / total

quartets_list = []
for quartet in quartets:
    quartets_list.append(f"{quartet} {quartets[quartet]}\n")

quartets_list.sort()

with open(output_quartets_file, "w") as fp:
    for quartet in quartets_list:
        fp.write(quartet)

print(f"Elapsed time: {(time.time()-start_time)} seconds")


        

