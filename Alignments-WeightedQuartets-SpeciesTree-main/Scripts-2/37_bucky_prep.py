import os, sys
import itertools

####################### COMMANDS ########################
# To perform bootstrap on 37.fasta, output written to RAxML_bootstrap.all
# raxmlHPC -m GTRGAMMA -n all -# 200 -s 37.fasta -p 12345 -b 12345


####################### GLOBAL VARIABLES ########################
configuration = '1X-200-1000'
# if (len(sys.argv) > 1):
#     configuration = sys.argv[1]

GENES = int(configuration.split('-')[1])
REPLICATES = 20

# base_folder = '/home/navidh86/Desktop/comparative_study/Alignments-WeightedQuartets-SpeciesTree-main'
base_folder = '/home/navid/comparative_study/Alignments-WeightedQuartets-SpeciesTree-main'
input_folder = f'{base_folder}/37-taxon/{configuration}'
output_folder = f'{base_folder}/37-taxon/{configuration}'
stree_output_folder = f"{base_folder}/Estimated-Species-Trees/37-taxon/{configuration}"
wqrts_output_folder = f'{base_folder}/37-taxon/Output/{configuration}'

os.makedirs(output_folder, exist_ok=True)
os.makedirs(stree_output_folder, exist_ok=True)
os.makedirs(wqrts_output_folder, exist_ok=True)

#------- default MrBayes parameters ---------#
nst   = 6;       # HKY model, nst=2. Use $nst = 6 for GTR
rates = "gamma"; # rate variation across sites: Gamma distribution
                  # use "invgamma" for Gamma+I
ngen=55000
# ngen = 110000
nruns=2
mbsumburnin=6 # should be 11
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
            gene_file = gene + (replicate_num-1) * GENES
            input_file_sequences_fasta = os.path.join(input_folder, f"R{replicate_num}/{gene_file}/{gene_file}.fasta")
            output_file_sequences_nexus = os.path.join(output_folder_rep, f"{gene_file}.nex")
            cmd = "python3 convert_formats.py -i {} -o {} -m1 fasta -m2 nexus".format(input_file_sequences_fasta, output_file_sequences_nexus)
            print(cmd)
            os.system(cmd)

def prepare_mr_bayes_commands():
    for replicate_num in range(1, REPLICATES+1):
    # for replicate_num in range(1, 1+1):
        output_folder_rep = os.path.join(output_folder, f"R{replicate_num}/all_nexus")
        os.makedirs(output_folder_rep, exist_ok=True)

        for gene in range(1, GENES+1):
        # for gene in range(1, 1+1):
            gene_file = gene + (replicate_num-1) * GENES
            nexus_sequence_file = os.path.join(output_folder_rep, f"{gene_file}.nex")
            nexus_command_file = os.path.join(output_folder_rep, f"mb{gene_file}.nex")
            mbOutFile = os.path.join(output_folder_rep, f"{gene_file}")

            # write the commands
            with open(nexus_command_file, "w") as fp:
                fp.write("#nexus\nbegin mrbayes;\nset autoclose=yes nowarn=yes;\n")
                fp.write(f"execute {gene_file}.nex;\n\n")
                # fp.write(f"execute {nexus_sequence_file};\n")
                fp.write(f"lset nst={nst} rates={rates};\n")
                fp.write("prset brlenspr=Unconstrained:Exp(50.0);\n")
                # above: prior mean of 1/50=0.02 for branch lengths.
                fp.write(f"mcmc nruns={nruns} temp={temp} ngen={ngen} burninfrac={mbburninfrac} Nchains={Nchains} ")
                fp.write(f"samplefreq={samplefreq} swapfreq={swapfreq} printfreq={printfreq} ")
                fp.write(f"mcmcdiagn=yes diagnfreq={diagnfreq} ")
                fp.write(f"filename={gene_file};\n")

                # consensus tree (this was added later)
                # "sumt nruns=$nruns filename=$mbOutFile burnin=$mbsumburnin;\n";
                fp.write(f"sumt nruns={nruns} filename={gene_file} burnin={mbsumburnin};\n")

                # fp.write(f"filename={mbOutFile};\n")
                fp.write("quit;\nend;\n")

def run_mr_bayes():
    curdir = os.getcwd()
    
    rep_start = int(sys.argv[1])
    rep_end = int(sys.argv[2])

    # for replicate_num in range(1, REPLICATES+1):
    for replicate_num in range(rep_start, rep_end+1):
        output_folder_rep = os.path.join(output_folder, f"R{replicate_num}/all_nexus")
        os.makedirs(output_folder_rep, exist_ok=True)
        os.chdir(output_folder_rep)

        # if replicate_num == 2:
        #     gene_start = 60
        # elif replicate_num == 12:
        #     gene_start = 51
        # else:
        #     gene_start = 1
        gene_start = 1

        for gene in range(gene_start, GENES+1):
        # for gene in range(gene_start, gene_start+1):
            gene_file = gene + (replicate_num-1) * GENES
            cmd = f"mb mb{gene_file}.nex >> {gene_file}.log"
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
            gene_file = gene + (replicate_num-1) * GENES
            cmd = f"mbsum -n {mbsumburnin} -o {gene_file}.in {gene_file}.run?.t"
            print(replicate_num, cmd)
            os.system(cmd)

            os.remove(f"{gene_file}.run1.t")
            os.remove(f"{gene_file}.run2.t")

            os.system(f"mv {gene_file}.in {output_folder_rep}")

    os.chdir(curdir)


taxa_list = ["GAL", "ORN", "PRO", "LOX", "ECH", "DAS", "CHO", "OTO", "MIC", "CAL", "PON", "GOR",
             "HOM", "PAN", "NEW", "TAR", "OCH", "ORY", "SPE", "DIP", "RAT", "MUS", "CAV", "TUP",
             "PTE", "MYO", "EQU", "FEL", "CAN", "SUS", "BOS", "TUR", "VIC", "ERI", "SOR", "MON",
             "MAC"]

def convert_to_numbers(tree):
    # need a fixed mapping for this, as all bootstrapped trees need to follow the same mapping
    for nmb in range(1, 38):
        tree = tree.replace(taxa_list[nmb-1], str(nmb))

    return tree

def prepare_bootstrap_input():
    for replicate_num in range(1, REPLICATES+1):
        print(f"Preparing replicate {replicate_num}")
        input_folder_rep = os.path.join(input_folder, f"R{replicate_num}")
        output_folder_rep = os.path.join(output_folder, f"R{replicate_num}/all_in_bootstrap")
        os.makedirs(output_folder_rep, exist_ok=True)

        for gene in range(1, GENES+1):
            print(f"Gene {gene}")
            input_file = os.path.join(input_folder_rep, 
                        f"{gene+(replicate_num-1)*GENES}/raxmlboot.gtrgamma/RAxML_bootstrap.all")
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
                for nmb in range(1, 37):
                    fp.write(f"       {nmb} {taxa_list[nmb-1]},\n")
                fp.write(f"       37 {taxa_list[36]};\n")

                for tree in counts.keys():
                    fp.write(f"{tree} {counts[tree]}\n")


alpha = '1-s500'
gt_type = 'all_in' # 'all_in_bootstrap'

def run_bucky():
    n = 100000 # number of MCMC updates
    c = 2 # number of chains
    k = 2 # number of independent runs

    curdir = os.getcwd()

    for replicate_num in range(1, REPLICATES+1):
        input_folder_rep = os.path.join(output_folder, f"R{replicate_num}/{gt_type}")
        output_file = os.path.join(input_folder_rep, f"a{alpha}")
        os.makedirs(input_folder_rep, exist_ok=True)
        os.chdir(input_folder_rep)

        if alpha == "inf":
            cmd = f"bucky -n {n} -k {k} -c {c} --use-independence-prior -o {output_file} {input_folder_rep}/*.in"
        else:
            cmd = f"bucky -n {n} -k {k} -c {c} -a 1 -o {output_file} {input_folder_rep}/*.in"
        
        print(replicate_num, cmd)
        os.system(cmd)

        os.system(f"rm {input_folder_rep}/*.cluster")
        os.system(f"rm {input_folder_rep}/*.gene")
        os.system(f"rm {input_folder_rep}/*.input")
        os.system(f"rm {input_folder_rep}/*.out")

    os.chdir(curdir)


mapping_int_str = {}
def translate():
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

        for num in range(37, 0, -1):
            population_tree = population_tree.replace(str(num), mapping_int_str[replicate_num][num])
            concordance_tree = concordance_tree.replace(str(num), mapping_int_str[replicate_num][num])

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
        input_folder = f"{base_folder}/37-taxon/{configuration}/R{replicate_num}/all_nexus"
        output_folder = f"{base_folder}/37-taxon/{configuration}/R{replicate_num}/mb-consensus-trees"
        
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

        for number in range(37, 0, -1):
            tree = tree.replace(str(number), mapping_int_str[replicate_num][number])

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

        for number in range(37, 0, -1):
            tree = tree.replace(str(number), mapping_int_str[replicate_num][number])

        with open(output_tree_file, "w") as fout:
            fout.write(tree)

        os.remove(temporary_tree_file)

# prepare_nexus_sequences()
# prepare_mr_bayes_commands()
# run_mr_bayes()
# move_consensus_trees()
# run_mbsum()
# prepare_bootstrap_input()

# # for g in ['all_in', 'all_in_bootstrap']:
for g in ['all_in']:
    # gt_type = g
    run_bucky()
    # translate()
    # extract_bucky_species_trees()
    # extract_bucky_species_quartets()
    
    # run_wqfm()
    # run_wqmc()
    pass

