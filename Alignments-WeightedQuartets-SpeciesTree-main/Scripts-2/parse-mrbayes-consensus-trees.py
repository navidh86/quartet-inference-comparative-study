import sys
from newick import loads

GENES = 200 if len(sys.argv) == 1 else int(sys.argv[1])

for gene in range(1, GENES+1):
    with open(f"{gene}.con.tre", "r") as fp:
        lines = fp.readlines()

    mapping = {}
    TRANSLATE = 0 # 0 -> waiting to translate, 1 -> translating; 2 -> looking for tree, 3 -> done
    for line in lines:
        # print(line)
        if TRANSLATE == 0:
            if "translate" in line:
                TRANSLATE = 1
        elif TRANSLATE == 1:
            if ";" in line:
                TRANSLATE = 2
            else:
                line = line.strip()
                tmp = line.split("\t")
                if tmp[1][-1] == ",":
                    tmp[1] = tmp[1][:-1]
                mapping[int(tmp[0])] = tmp[1]
        elif TRANSLATE == 2:
            tree = line[line.find("("):line.find(";")+1]
            break

    print(mapping)
        
    print("################ BEFORE ###################")
    print(tree)

    old_tree = tree + ""

    # remove branch lengths
    while ":" in tree:
        part1 = tree[:tree.find(":")]
        part2 = tree[tree.find("[&length")+8:]
        tree = f"{part1}[{part2}"


    # remove support
    while "[" in tree:
        part1 = tree[:tree.find("[")]
        part2 = tree[tree.find("]")+1:]
        tree = f"{part1}{part2}"


    print("\n\n\n\n################ AFTER ###################")
    print(tree)
    taxa = [int(taxon) for taxon in tree.replace("(", "").replace(")", "").replace(";", "").split(",")]

    # translate
    for label in range(37, 0, -1):
        tree = tree.replace(str(label), mapping[label])

    print("\n\n\n\n################ FINALLY ###################")
    print(tree)

    print("\n\n\n\n################ Draw Probability (Support) ###################")
    supports = []
    kwl = "prob(percent)"
    kwr = "prob+-sd"
    while kwl in old_tree:
        support = old_tree[old_tree.find(kwl)+len(kwl)+2:old_tree.find(kwr)-2]
        supports.append(float(support))
        old_tree = old_tree.replace(kwl, "done", 1)
        old_tree = old_tree.replace(kwr, "done", 1)

    print(taxa)
    print(supports)

    cnt = 0; idx = 0
    while idx < len(tree):
        t = tree[idx]
        if t == ",":
            cnt += 1
        elif t == ")":
            cnt += 1
            if cnt < len(supports):
                tree = f"{tree[:idx+1]}{supports[cnt]}{tree[idx+1:]}"

        idx += 1

    print(tree)

    with open(f"{gene}.mb.tre", "w") as fp:
        fp.write(tree + "\n")
