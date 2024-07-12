import os, sys
import itertools

configuration = 'lower-ILS'
if len(sys.argv) > 1:
    configuration = sys.argv[1]
REPLICATES = 20
GENES = 200
base_folder = '/home/navid/comparative_study/Alignments-WeightedQuartets-SpeciesTree-main'
input_folder = f'{base_folder}/10-taxon/{configuration}/estimated-genetrees'
output_folder = f'{base_folder}/10-taxon/{configuration}/estimated-genetrees'
stree_output_folder = f"{base_folder}/Estimated-Species-Trees/10-taxon/{configuration}"
wqrts_output_folder = f'{base_folder}/10-taxon/Output/{configuration}'

#------- default MrBayes parameters ---------#
nst   = 6;       # HKY model, nst=2. Use $nst = 6 for GTR
rates = "gamma"; # rate variation across sites: Gamma distribution
                  # use "invgamma" for Gamma+I
ngen=55000
nruns=2
mbsumburnin=6
mbburninfrac=5000/55000  #100/1100
diagnfreq=10000
# samplefreq=1000
samplefreq=500
printfreq=100
temp=0.04
Nchains=4
swapfreq=10
# Nswaps=1

# convert the gene sequences to nexus and move them to the proper folder
def prepare_nexus_sequences():
    for replicate_num in range(1, REPLICATES+1):
        output_folder_rep = os.path.join(output_folder, f"R{replicate_num}/all_nexus")
        os.makedirs(output_folder_rep, exist_ok=True)

        for gene in range(1, GENES+1):
            input_file_sequences_fasta = os.path.join(input_folder, f"R{replicate_num}/{gene}/truegene.fasta")
            output_file_sequences_nexus = os.path.join(output_folder_rep, f"{gene}.nex")
            cmd = "python3 convert_formats.py -i {} -o {} -m1 fasta -m2 nexus".format(input_file_sequences_fasta, output_file_sequences_nexus)
            print(cmd)
            os.system(cmd)

def prepare_mr_bayes_commands():
    for replicate_num in range(1, REPLICATES+1):
        output_folder_rep = os.path.join(output_folder, f"R{replicate_num}/all_nexus")
        os.makedirs(output_folder_rep, exist_ok=True)

        for gene in range(1, GENES+1):
            nexus_sequence_file = os.path.join(output_folder_rep, f"{gene}.nex")
            nexus_command_file = os.path.join(output_folder_rep, f"mb{gene}.nex")
            mbOutFile = os.path.join(output_folder_rep, f"{gene}")

            # write the commands
            with open(nexus_command_file, "w") as fp:
                fp.write("#nexus\nbegin mrbayes;\nset autoclose=yes nowarn=yes;\n")
                fp.write(f"execute {gene}.nex;\n\n")
                # fp.write(f"execute {nexus_sequence_file};\n")
                fp.write(f"lset nst={nst} rates={rates};\n")
                fp.write("prset brlenspr=Unconstrained:Exp(50.0);\n")
                # above: prior mean of 1/50=0.02 for branch lengths.
                fp.write(f"mcmc nruns={nruns} temp={temp} ngen={ngen} burninfrac={mbburninfrac} Nchains={Nchains} ")
                fp.write(f"samplefreq={samplefreq} swapfreq={swapfreq} printfreq={printfreq} ")
                fp.write(f"mcmcdiagn=yes diagnfreq={diagnfreq} ")
                fp.write(f"filename={gene};\n")
                fp.write(f"sumt nruns={nruns} filename={gene} burnin={mbsumburnin};\n")
                # fp.write(f"filename={mbOutFile};\n")
                fp.write("quit;\nend;\n")

def run_mr_bayes():
    curdir = os.getcwd()
    
    for replicate_num in range(1, REPLICATES+1):
        output_folder_rep = os.path.join(output_folder, f"R{replicate_num}/all_nexus")
        os.makedirs(output_folder_rep, exist_ok=True)
        os.chdir(output_folder_rep)

        # if replicate_num == 1:
        #     gene_start = 232
        # else:
        #     gene_start = 1

        for gene in range(1, GENES+1):
            cmd = f"mb mb{gene}.nex >> {gene}.log"
            print(replicate_num, cmd)
            os.system(cmd)

            os.remove(f"{gene_file}.run1.p")
            os.remove(f"{gene_file}.run2.p")
            os.remove(f"{gene_file}.mcmc")
            os.remove(f"{gene_file}.ckp")
            os.remove(f"{gene_file}.log")
            os.remove(f"{gene_file}.parts")
            os.remove(f"{gene_file}.trprobs")
            os.remove(f"{gene_file}.tstat")
            os.remove(f"{gene_file}.vstat")

    os.chdir(curdir)

def run_mbsum():
    curdir = os.getcwd()
    
    for replicate_num in range(1, REPLICATES+1):
        input_folder_rep = os.path.join(output_folder, f"R{replicate_num}/all_nexus")
        output_folder_rep = os.path.join(output_folder, f"R{replicate_num}/all_in")
        os.makedirs(input_folder_rep, exist_ok=True)
        os.makedirs(output_folder_rep, exist_ok=True)
        os.chdir(input_folder_rep)

        for gene in range(1, GENES+1):
            cmd = f"mbsum -n {mbsumburnin} -o {gene}.in {gene}.run?.t"
            print(replicate_num, cmd)
            os.system(cmd)

            os.remove(f"{gene}.run1.t")
            os.remove(f"{gene}.run2.t")

            os.system(f"mv {gene}.in {output_folder_rep}")

    os.chdir(curdir)

def convert_to_numbers(tree):
    # need to change in two levels
    # first convert 0-10 to A-K
    for nmb in range(10, -1, -1):
        tree = tree.replace(str(nmb), chr(ord('A')+nmb))
    
    # now replace the characters A-K by 1-11
    for nmb in range(0, 11):
        tree = tree.replace(chr(ord('A')+nmb), str(nmb+1))

    return tree

def prepare_bootstrap_input():
    for replicate_num in range(1, REPLICATES+1):
        input_folder_rep = os.path.join(output_folder, f"R{replicate_num}")
        output_folder_rep = os.path.join(output_folder, f"R{replicate_num}/all_in_bootstrap")
        os.makedirs(output_folder_rep, exist_ok=True)

        for gene in range(1, GENES+1):
            input_file = os.path.join(input_folder_rep, f"{gene}/RAxML_bootstrap.all")
            output_file = os.path.join(output_folder_rep, f"{gene}.in")

            counts = {}
            trees = []
            with open(input_file, "r") as fp:
                lines = fp.readlines()

            for line in lines:
                tree = convert_to_numbers(line[:-1])
                if (tree in counts):
                    counts[tree] += 1
                else:
                    counts[tree] = 1

            with open(output_file, "w") as fp:
                fp.write("translate\n")
                for c in range(0, 10):
                    fp.write(f"       {c+1} {c},\n")
                fp.write("       11 10;\n")

                for tree in counts.keys():
                    fp.write(f"{tree} {counts[tree]}\n")

def rename_files():
    methods = ['bucky-conc-mrbayes.tre',
               'bucky-conc-RAxML.tre',
               'bucky-mrbayes.tre',
               'bucky-RAxML.tre',
               'wQFM-bucky-mrbayes.tre',
               'wQFM-bucky-RAxML.tre',
               'wQMC-bucky-mrbayes.tre',
               'wQMC-bucky-RAxML.tre'
               ]
    for replicate_num in range(1, REPLICATES+1):
        for method in methods:
            src_name = f"{stree_output_folder}/R{replicate_num}-a1-{method}"
            dest_name = f"{stree_output_folder}/R{replicate_num}-ainf-{method}"

            os.rename(src_name, dest_name)

if len(sys.argv) > 2:
    alpha = sys.argv[2]
else:
    alpha = '1-s500'

gt_type = 'all_in' # 'all_in_bootstrap'

def run_bucky():
    n = 100000 # number of MCMC updates
    c = 2 # number of chains
    k = 2 # number of independent runs

    curdir = os.getcwd()

    # for replicate_num in range(1, REPLICATES+1):
    for replicate_num in range(1, 1+1):
        input_folder_rep = os.path.join(output_folder, f"R{replicate_num}/{gt_type}")
        output_file = os.path.join(input_folder_rep, f"a{alpha}-test")
        os.makedirs(input_folder_rep, exist_ok=True)
        os.chdir(input_folder_rep)

        if alpha == "inf":
            cmd = f"bucky -n {n} -k {k} -c {c} --use-independence-prior -o {output_file} {input_folder_rep}/*.in"
        else:
            # cmd = f"bucky -n {n} -k {k} -c {c} -a 1 -o {output_file} {input_folder_rep}/*.in"
            cmd = f"bucky -n {n} -k {k} -c {c} -a 1 --create-single-file -o {output_file} {input_folder_rep}/*.in"
        
        print(replicate_num, cmd)
        os.system(cmd)

        # os.system(f"rm {input_folder_rep}/*.cluster")
        # os.system(f"rm {input_folder_rep}/*.gene")
        # os.system(f"rm {input_folder_rep}/*.input")
        # os.system(f"rm {input_folder_rep}/*.out")

    os.chdir(curdir)


mapping_int_str = {}
def translate():
    # reset
    global mapping_int_str
    mapping_int_str = {}

    for replicate_num in range(1, REPLICATES+1):
        concordance_file = os.path.join(output_folder, f"R{replicate_num}/{gt_type}/a{alpha}.concordance")
        # Change back to original names
        # create a integer to string mapping first
        with open(concordance_file, "r") as fin:
            lines = [l.strip() for l in fin.readlines()]

        def get_taxa_int_str(line):
            line = line.replace(",", "") # remove COMMAs
            line = line.replace(";", "") # remove COMMAs
            arr = line.split()
            return int(arr[0]), str(arr[1])
        
        mapping_int_str[replicate_num] = {}

        state_current = "UNINITIATED"
        for line in lines:
            if state_current == "TRANSLATE" and line != "":
                taxa_int, taxa_str = get_taxa_int_str(line)
                mapping_int_str[replicate_num][taxa_int] = taxa_str
            elif state_current == "UNINITIATED" and line == "translate":
                state_current = "TRANSLATE"
            elif state_current == "TRANSLATE" and line == "":
                state_current = "END"
                break

def custom_replace(tree, mapping):
    # can't use the mapping directly, as both key and value have common numbers
    # first, we create two maps, one maps keys to A-K, the other maps A-K to the values
    map1 = {}; map2 = {}
    ch = ord('A')
    for key in mapping:
        map1[key] = chr(ch)
        map2[chr(ch)] = mapping[key]
        ch += 1

    # now replace the keys (1-11) with map1
    for nmb in range(11, 0, -1):
        tree = tree.replace(str(nmb), map1[nmb])

    # now replace them using map2
    for nmb in range(0, 11):
        tree = tree.replace(chr(ord('A')+nmb), map2[chr(ord('A')+nmb)])

    return tree

def extract_bucky_species_trees():
    for replicate_num in range(1, REPLICATES+1):
        concordance_file = os.path.join(output_folder, f"R{replicate_num}/{gt_type}/a{alpha}.concordance")
        if gt_type == 'all_in_bootstrap':
            output_population_tree_file = os.path.join(stree_output_folder, f"R{replicate_num}-a{alpha}-bucky-RAxML.tre")
            output_concordance_tree_file = os.path.join(stree_output_folder, f"R{replicate_num}-a{alpha}-bucky-conc-RAxML.tre")
        else:
            output_population_tree_file = os.path.join(stree_output_folder, f"R{replicate_num}-a{alpha}-bucky-mrbayes.tre")
            output_concordance_tree_file = os.path.join(stree_output_folder, f"R{replicate_num}-a{alpha}-bucky-conc-mrbayes.tre")

        # Change back to original names
        # create a integer to string mapping first
        with open(concordance_file, "r") as fin:
            lines = [l.strip() for l in fin.readlines()]

        def get_taxa_int_str(line):
            line = line.replace(",", "") # remove COMMAs
            line = line.replace(";", "") # remove COMMAs
            arr = line.split()
            return int(arr[0]), str(arr[1])
        
        # mapping_int_str = {} # using the global one now

        state_current = "TREESEEK"
        for line in lines:
            if state_current == "TREESEEK" and line == "Population Tree:":
                state_current = "TREEREAD"
            elif state_current == "TREEREAD":
                population_tree = line
                state_current = "TREESEEKCONC"
            elif state_current == "TREESEEKCONC" and line == "Primary Concordance Tree Topology:":
                state_current = "TREEREADCONC"
            elif state_current == "TREEREADCONC":
                concordance_tree = line
                break

        # for num in range(11, 0, -1):
        #     population_tree = population_tree.replace(str(num), mapping_int_str[replicate_num][num])
        #     concordance_tree = concordance_tree.replace(str(num), mapping_int_str[replicate_num][num])

        population_tree = custom_replace(population_tree, mapping_int_str[replicate_num])
        concordance_tree = custom_replace(concordance_tree, mapping_int_str[replicate_num])

        with open(output_population_tree_file, "w") as fout:
            fout.write(population_tree)

        with open(output_concordance_tree_file, "w") as fout:
            fout.write(concordance_tree)
                
def extract_bucky_species_quartets():
    # Four-way partitions
    for replicate_num in range(1, REPLICATES+1):
        # concordance_file = os.path.join(output_folder, f"R{replicate_num}/all_in/a{alpha}.concordance")
        concordance_file = os.path.join(output_folder, f"R{replicate_num}/{gt_type}/a{alpha}.concordance")
        if gt_type == 'all_in_bootstrap':
            output_wqrts_file_wqfm = os.path.join(wqrts_output_folder, f"R{replicate_num}-a{alpha}-bucky-RAxML.wqrts")
            output_wqrts_file_wqmc = os.path.join(wqrts_output_folder, f"R{replicate_num}-a{alpha}-bucky-RAxML.wqmc.wqrts")
        else:
            output_wqrts_file_wqfm = os.path.join(wqrts_output_folder, f"R{replicate_num}-a{alpha}-bucky-mrbayes.wqrts")
            output_wqrts_file_wqmc = os.path.join(wqrts_output_folder, f"R{replicate_num}-a{alpha}-bucky-mrbayes.wqmc.wqrts")
     
        with open(concordance_file, "r") as fin:
            lines = [l.strip() for l in fin.readlines()]

        quartets = []
        cfs = []
        cus = []

        state_current = "QSEEK"
        for line in lines:
            if state_current == "QSEEK" and "Four-way partitions" in line:
                state_current = "QREAD"
            elif state_current == "QREAD" and line == "":
                state_current = "END"
                break
            elif state_current == "QREAD":
                # {1,3,5,6,7,8,9,10,11,12,13,14; 4|2; 15}	0.388, 0.086,
                quartets.append(line.split("}")[0][1:])
                temp = line.split("}\t")[1].split(" ")
                cfs.append(temp[0][:-1])
                cus.append(temp[1][:-1])

        # 1,3,5,6,7,8,9,10,11,12,13,14; 4|2; 15
        wqrts = []; wqmc = []
        # ((A,B),(C,D));
        for idx, quartet in enumerate(quartets):
            t = [[] for x in range(4)]
            temp = quartet.split("|")
            temp1 = temp[0].split("; ")
            t[0] = temp1[0].split(",")
            t[1] = temp1[1].split(",")
            
            temp2 = temp[1].split("; ")
            t[2] = temp2[0].split(",")
            t[3] = temp2[1].split(",")

            # print(t)
            # wqs = list(itertools.product(*t))
            # print(wqs)
            # # for wq in :
            #     # print(wq)
            # break


            # ((A,B),(C,D)); 41
            for t0 in t[0]:
                for t1 in t[1]:
                    for t2 in t[2]:
                        for t3 in t[3]:
                            qrt = f"(({t0},{t1}),({t2},{t3})); {cfs[idx]}"
                            wqmc_qrt = f"{t0},{t1}|{t2},{t3}:{cfs[idx]}"
                            wqrts.append(qrt)
                            wqmc.append(wqmc_qrt)

        with open(output_wqrts_file_wqfm, "w") as fp:
            for wqrt in wqrts:
                fp.write(wqrt + "\n")

        with open(output_wqrts_file_wqmc, "w") as fp:
            for wqrt in wqmc:
                fp.write(wqrt + "\n")

def move_consensus_trees():
    # For wastral input
    for replicate_num in range(1, REPLICATES+1):
        print(f"Replicate {replicate_num}......")
        input_folder = f"{base_folder}/10-taxon/{configuration}/estimated-genetrees/R{replicate_num}/all_nexus"
        output_folder = f"{base_folder}/10-taxon/{configuration}/estimated-genetrees/R{replicate_num}/mb-consensus-trees"
        
        os.makedirs(output_folder, exist_ok=True)

        cmd = f"mv {input_folder}/*.con.tre {output_folder}"
        print(cmd)
        os.system(cmd) 

def run_wqfm():
    for replicate_num in range(1, REPLICATES+1):
        print(f"Replicate {replicate_num}......")
        # input_wqrts_file = os.path.join(wqrts_output_folder, f"R{replicate_num}-bucky.wqrts")
        if gt_type == 'all_in_bootstrap':
            input_wqrts_file = os.path.join(wqrts_output_folder, f"R{replicate_num}-a{alpha}-bucky-RAxML.wqrts")
            wQFM_tree_file = os.path.join(stree_output_folder, f"R{replicate_num}-a{alpha}-wQFM-bucky-RAxML.tre")
        else:
            input_wqrts_file = os.path.join(wqrts_output_folder, f"R{replicate_num}-a{alpha}-bucky-mrbayes.wqrts")
            wQFM_tree_file = os.path.join(stree_output_folder, f"R{replicate_num}-a{alpha}-wQFM-bucky-mrbayes.tre")

        temporary_tree_file = os.path.join(stree_output_folder, f"R{replicate_num}-wQFM-bucky.temp")

        cmd = f"java -jar wQFM-v1.4.jar -i {input_wqrts_file} -o {temporary_tree_file}"
        print(cmd)
        os.system(cmd)

        # replace numbers with taxa
        with open(temporary_tree_file, "r") as fin:
            tree = fin.readline()

        # for number in range(11, 0, -1):
        #     tree = tree.replace(str(number), mapping_int_str[replicate_num][number])
        tree = custom_replace(tree, mapping_int_str[replicate_num])

        with open(wQFM_tree_file, "w") as fout:
            fout.write(tree)

        os.remove(temporary_tree_file)
                

def run_wqmc():
     for replicate_num in range(1, REPLICATES+1):
        print(f"Replicate {replicate_num}......")
        if gt_type == 'all_in_bootstrap':
            input_wqrts_file = os.path.join(wqrts_output_folder, f"R{replicate_num}-a{alpha}-bucky-RAxML.wqmc.wqrts")
            output_tree_file = os.path.join(stree_output_folder, f"R{replicate_num}-a{alpha}-wQMC-bucky-RAxML.tre")
        else:
            input_wqrts_file = os.path.join(wqrts_output_folder, f"R{replicate_num}-a{alpha}-bucky-mrbayes.wqmc.wqrts")
            output_tree_file = os.path.join(stree_output_folder, f"R{replicate_num}-a{alpha}-wQMC-bucky-mrbayes.tre")

        temporary_tree_file = os.path.join(stree_output_folder, f"R{replicate_num}-wQMC-bucky.temp")
        cmd = "./max-cut-tree qrtt={} weights=on otre={}".format(input_wqrts_file, temporary_tree_file)
        print(cmd)
        os.system(cmd)

		# replace numbers with taxa
        with open(temporary_tree_file, "r") as fin:
            tree = fin.readline()

        # for number in range(11, 0, -1):
        #     tree = tree.replace(str(number), mapping_int_str[replicate_num][number])
        tree = custom_replace(tree, mapping_int_str[replicate_num])

        with open(output_tree_file, "w") as fout:
            fout.write(tree)

        os.remove(temporary_tree_file)

if __name__ == "__main__":
    # prepare_nexus_sequences()
    # prepare_mr_bayes_commands()
    # run_mr_bayes()
    # move_consensus_trees()
    # run_mbsum()
    # prepare_bootstrap_input()

    # for g in ['all_in', 'all_in_bootstrap']:
    for g in ['all_in']:
        gt_type = g
        run_bucky()
        translate()
        # extract_bucky_species_trees()
        # extract_bucky_species_quartets()
        
        # run_wqfm()
        # run_wqmc()

    # rename_files()
