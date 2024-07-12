import os
import sys

####################### GLOBAL VARIABLES ########################
# configurations = ['1X-200-250', '0.5X-200-500', '1X-800-1000']
# cidx = 1 ##### CHANGE THIS ONLY
# configuration = configurations[cidx]
# configuration ='1X-200-500'
if len(sys.argv) > 1:
	configuration = sys.argv[1]
else:
	print("Please enter config name")
	exit(-1)

GENES = int(configuration.split("-")[1])
REPLICATES = 20

####################### FOLDERS ########################
# print(os.path.dirname(os.path.dirname((os.path.abspath(__file__)))))
# exit(0)
base_folder = os.path.dirname(os.path.dirname((os.path.abspath(__file__))))
base_folder = '/home/navid/comparative_study/Alignments-WeightedQuartets-SpeciesTree-main'
input_folder = f'{base_folder}/37-taxon/{configuration}'
output_folder = f'{base_folder}/37-taxon/Output/{configuration}'
astral_path = f"{base_folder}/Scripts-2/Astral.5.7.8/Astral/astral.5.7.8.jar"
stree_output_folder = f"{base_folder}/Estimated-Species-Trees/37-taxon/{configuration}"
fp_fn_folder = f"{base_folder}/Estimated-Species-Trees/FP_FN/37-taxon/{configuration}"
rf_folder = f"{base_folder}/Estimated-Species-Trees/RF-rates/37-taxon/{configuration}"

# create the folders if not available
os.makedirs(input_folder, exist_ok=True)
os.makedirs(output_folder, exist_ok=True)
os.makedirs(stree_output_folder, exist_ok=True)
os.makedirs(fp_fn_folder, exist_ok=True)
os.makedirs(rf_folder, exist_ok=True)

status_file = "status_37.txt"

def clear_status_file():
	with open(status_file, 'w'):
		pass

def form_fasta_text_files():
	for replicate_num in range(1, REPLICATES+1):
		replicate_folder = os.path.join(input_folder, f"R{replicate_num}")
		text_file_path = os.path.join(output_folder, f"R{replicate_num}-Fasta-Files.txt")

		print(f"Running for {replicate_folder}, text_file_path = {text_file_path}")
		
		with open(text_file_path, mode='w') as fout:
			for folder_itr in range(1, GENES+1):
				folder = "{}".format(folder_itr + (replicate_num-1)*GENES)
				fasta_file = os.path.join(replicate_folder, folder, f"{folder}.fasta")

				fout.write(fasta_file)
				fout.write("\n")

	with open(status_file, "a") as fout:
		fout.write("fasta files formed\n")

def concat_sequences():
	for replicate_num in range(1, REPLICATES+1):
		print("Running for replicate", replicate_num)
		text_file_path = os.path.join(output_folder, f"R{replicate_num}-Fasta-Files.txt")
		output_file_sequences_fasta = os.path.join(output_folder, f"R{replicate_num}-Sequences-Concatenated.fasta")
		output_file_sequences_nexus = os.path.join(output_folder, f"R{replicate_num}-Sequences-Concatenated.nex")
		
		###### CONCATENATE SEQUENCES
		cmd = "python3 concatenate_sequences.py -i {} -o {}".format(text_file_path, output_file_sequences_fasta)  
		## -n NOT_GIVEN_FOR_NOW
		print(cmd)
		os.system(cmd)

		##### CONVERT TO NEXUS FORMAT
		cmd = "python3 convert_formats.py -i {} -o {} -m1 fasta -m2 nexus".format(output_file_sequences_fasta, output_file_sequences_nexus)
		print(cmd)
		os.system(cmd)

	with open(status_file, "a") as fout:
		fout.write("sequences conacatenated\n")

def run_svd_wqfm():
	for mode in ["SVD-reciprocal", "SVD-exponential"]:
		for replicate_num in range(1, REPLICATES+1):
			input_file_sequences_nexus = os.path.join(output_folder, f"R{replicate_num}-Sequences-Concatenated.nex")
			output_file_wqrts = os.path.join(output_folder, f"R{replicate_num}-{mode}.wqrts")
			cmd = "python3 generate_SVD_wqrts.py -i {} -o {} -m {}".format(input_file_sequences_nexus, output_file_wqrts, mode[4:])

			print(cmd)
			os.system(cmd)

			# RUN wQFM
			wQFM_file = os.path.join(stree_output_folder, f"R{replicate_num}-wQFM-{mode}.tre")
			cmd = "java -jar wQFM-v1.4.jar -i {} -o {}".format(output_file_wqrts, wQFM_file)

			print(cmd)
			os.system(cmd)

	with open(status_file, "a") as fout:
		fout.write("SVD wQFM done\n")

def run_SVD_output_tree():
	for replicate_num in range(1, REPLICATES+1):
		sequence_file = os.path.join(output_folder, f"R{replicate_num}-Sequences-Concatenated.nex")
		output_tree_file = os.path.join(stree_output_folder, f"R{replicate_num}-SVD.tre")

		cmd = "python3 generate_SVD_tree.py -i {} -o {}".format(sequence_file, output_tree_file)
		print(cmd)
		os.system(cmd)

	with open(status_file, "a") as fout:
		fout.write("SVD output done\n")


def run_svd_wqmc():
	for mode in ["SVD-reciprocal", "SVD-exponential"]:
		for replicate_num in range(1, REPLICATES+1):
			input_file_sequences_nexus = os.path.join(output_folder, f"R{replicate_num}-Sequences-Concatenated.nex")
			output_file_wqrts = os.path.join(output_folder, f"R{replicate_num}-{mode}-wqmc.wqrts")
			taxa_mapping_file = os.path.join(output_folder, f"R{replicate_num}-{mode}-wqmc.mapping")

			cmd = "python3 generate_SVD_wqrts_wqmc.py -i {} -o {} -t {} -m {}".format(input_file_sequences_nexus, output_file_wqrts, taxa_mapping_file, mode[4:])

			print(cmd)
			os.system(cmd)

			# RUN wQMC
			output_tree_file = os.path.join(stree_output_folder, f"R{replicate_num}-wQMC-{mode}.tre")
			temporary_tree_file = os.path.join(stree_output_folder, f"R{replicate_num}-wQMC-{mode}.temp")
			cmd = "./max-cut-tree qrtt={} weights=on otre={}".format(output_file_wqrts, temporary_tree_file)
			print(cmd)
			os.system(cmd)

			# Change back to original names
			# create a integer to string mapping first
			with open(taxa_mapping_file, mode='r') as fin:
				lines = [l.strip() for l in fin.readlines()]

			##### states = ["UNINITIATED", "INITIATED", "TERMINATED"]
			state_current = "UNINITIATED"

			def get_taxa_int_str(line):
				line = line.replace(",", "") # remove COMMAs
				arr = line.split()
				return int(arr[0]), str(arr[1])

			mapping_int_str = {}

			for line in lines:
				if state_current == "INITIATED" and line != ";":
					taxa_int, taxa_str = get_taxa_int_str(line)
					mapping_int_str[taxa_int] = taxa_str

				if line == "Translate":
					state_current = "INITIATED" ## Enter DFA
				elif line == ";":
					state_current = "TERMINATED" ## Exit DFA


			# Replace
			with open(temporary_tree_file, "r") as fin:
				tree = fin.readline()

			for num in range(37, 0, -1):
				tree = tree.replace(str(num), mapping_int_str[num])

			with open(output_tree_file, "w") as fout:
				fout.write(tree)

			os.remove(temporary_tree_file)
			os.remove(taxa_mapping_file)

	with open(status_file, "a") as fout:
		fout.write("SVD-wQMC done\n")

def run_RaXML_wqrts_wQFM():
	for mode in ["RAxML"]:
		for replicate_num in range(1, REPLICATES):
			sequence_file = os.path.join(output_folder, f"R{replicate_num}-Sequences-Concatenated.nex")
			output_wqrt_file = os.path.join(output_folder, f"R{replicate_num}-{mode}.wqrts")
			wQFM_file = os.path.join(stree_output_folder, f"R{replicate_num}-wQFM-RAxML.tre")
						
			## Run RAxML
			cmd = "python3 raxml_runner.py {} {}".format(sequence_file, output_wqrt_file)
			print(cmd)
			os.system(cmd)

			## Generate wQFM
			cmd = "java -jar wQFM-v1.4.jar -i {} -o {}".format(output_wqrt_file, wQFM_file)
			print(cmd)
			os.system(cmd)
			os.remove("RAxML_info.TEST")

def concat_gene_trees():
	for replicate_num in range(1, REPLICATES+1):
		print(f"Replicate {replicate_num}......")
		combined_gene_tree_file = os.path.join(input_folder, f'R{replicate_num}/all_gt.tre')
		with open(combined_gene_tree_file, "w") as fout:
			for gene in range(1, GENES+1):
				gene_tree_file = os.path.join(input_folder, f'R{replicate_num}/{gene+GENES*(replicate_num-1)}/raxmlboot.gtrgamma/RAxML_bipartitions.final.f200')

				with open(gene_tree_file, "r") as fin:
					tree = fin.readline()

				fout.write(tree)
	with open(status_file, "a") as fout:
		fout.write("concat gene trees done\n")

def concat_bootstrap_gene_trees():
	# this will concat all 200 bootstrapped gene trees for each replicate
	for replicate_num in range(1, REPLICATES+1):
		print(f"Replicate {replicate_num}......")
		combined_gene_tree_file = os.path.join(input_folder, f'R{replicate_num}/all_boot_gt.tre')
		if os.path.exists(combined_gene_tree_file):
			os.remove(combined_gene_tree_file)
		for gene in range(1, GENES+1):
			gene_tree_file = os.path.join(input_folder, f'R{replicate_num}/{gene+GENES*(replicate_num-1)}/raxmlboot.gtrgamma/RAxML_bootstrap.all')

			cmd = f"cat {gene_tree_file} >> {combined_gene_tree_file}"
			os.system(cmd)

	with open(status_file, "a") as fout:
		fout.write("concat boot gene trees done\n")

def concat_bucky_boot_gene_trees():
	for replicate_num in range(1, REPLICATES+1):
		print(f"Replicate {replicate_num}......")
		combined_gene_tree_file = os.path.join(input_folder, f'R{replicate_num}/all_bucky_boot.tre')
		if os.path.exists(combined_gene_tree_file):
			os.remove(combined_gene_tree_file)
		for gene in range(1, GENES+1):
			g = (replicate_num-1)*GENES + gene
			gene_tree_file = os.path.join(input_folder, f'R{replicate_num}/all_in/{g}.in')

			with open(gene_tree_file, "r") as fin:
				lines = [l.strip() for l in fin.readlines()]

			def get_taxa_int_str(line):
				line = line.replace(",", "") # remove COMMAs
				line = line.replace(";", "") # remove COMMAs
				arr = line.split()
				return int(arr[0]), str(arr[1])
			
			# read the mapping
			mapping_int_str = {}
			for line in lines[1:38]:
				taxa_int, taxa_str = get_taxa_int_str(line)
				mapping_int_str[taxa_int] = taxa_str
			
			# write all trees (remapped) to the output file
			with open(combined_gene_tree_file, "a") as fp:
				for line in lines[38:]:
					tree, count = line.split()
					# remap
					for nmb in range(37, 0, -1):
						tree = tree.replace(str(nmb), mapping_int_str[nmb])
					[fp.write(tree + "\n") for x in range(int(count))]

	with open(status_file, "a") as fout:
		fout.write("concat bucky boot gene trees done\n")

def generate_gtf_wqrts(boot=False, bucky_boot=False, del_boot=False):
	for replicate_num in range(16, 20+1):
	# for replicate_num in range(1, REPLICATES+1):
		print(f"Replicate {replicate_num}......")
		if boot:
			combined_gene_tree_file = os.path.join(input_folder, f'R{replicate_num}/all_boot_gt.tre')
			output_wqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF-boot.wqrts")
		elif bucky_boot:
			combined_gene_tree_file = os.path.join(input_folder, f'R{replicate_num}/all_bucky_boot.tre')
			output_wqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF-bucky-boot.wqrts")
		else:
			combined_gene_tree_file = os.path.join(input_folder, f'R{replicate_num}/all_gt.tre')
			output_wqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF.wqrts")

		cmd = f"python3 generate_wqrts.py {combined_gene_tree_file} {output_wqrts_file}"
		print(cmd)
		os.system(cmd)

		if boot and del_boot:
			os.remove(combined_gene_tree_file)

	with open(status_file, "a") as fout:
		fout.write("generating gtf qrts done\n")

def generate_distance_based_gtf_wqrts(boot=False, del_boot=False):
	for replicate_num in range(1, REPLICATES+1):
	# for replicate_num in range(1, 5+1):
		print(f"Replicate {replicate_num}......")
		if boot:
			combined_gene_tree_file = os.path.join(input_folder, f'R{replicate_num}/all_boot_gt.tre')
			output_wqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF-boot.dwqrts")
		else:
			combined_gene_tree_file = os.path.join(input_folder, f'R{replicate_num}/all_gt.tre')
			output_wqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF.dwqrts")

		cmd = f"python3 generate_distance_based_quartets.py {combined_gene_tree_file} {output_wqrts_file}"
		print(cmd)
		os.system(cmd)

		if boot and del_boot:
			os.remove(combined_gene_tree_file)

	with open(status_file, "a") as fout:
		fout.write("generating distance based gtf qrts done\n")

def generate_dominant_quartets():
	for replicate_num in range(1, REPLICATES+1):
		print(f"Replicate {replicate_num}......")
		
		input_wqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF.wqrts")
		output_dqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF.dqrts")
		output_udqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF.udqrts")

		cmd = f"java -jar generateBestWQrts.jar {input_wqrts_file} {output_dqrts_file}"
		print(cmd)
		os.system(cmd)
		cmd = f"awk '{{print $1\" \"1}}' {output_dqrts_file} > {output_udqrts_file}"
		print(cmd)
		os.system(cmd)

	with open(status_file, "a") as fout:
		fout.write("generating dominant qrts done\n")

def adjust_quartets(w1=0.8, w2=0.2, boot=False, del_boot=False):
	for replicate_num in range(1, REPLICATES+1):
		print(f"Replicate {replicate_num}......")
		if boot:
			input_wqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF-boot.wqrts")
			output_wqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF-boot.awqrts")
		else:
			input_wqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF.wqrts")
			output_wqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF.awqrts")

		cmd = f"python3 quartet_correction/adjust_weighted.py {input_wqrts_file} {output_wqrts_file} {w1} {w2}"
		print(cmd)
		os.system(cmd)


	with open(status_file, "a") as fout:
		fout.write("correction of qrts done\n")

def exec_wqfm(quartet_file, tree_file):
	cmd = f"java -jar wQFM-v1.4.jar -i {quartet_file} -o {tree_file}"
	print(cmd)
	os.system(cmd)

def run_wqfm_gtf(adjusted=False, boot=False, distance=False, dominant=False, bucky_boot=False):
	for replicate_num in range(1, REPLICATES+1):
		print(f"Replicate {replicate_num}......")
		if dominant:
			input_dqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF.dqrts")
			wQFM_file = os.path.join(stree_output_folder, f"R{replicate_num}-wQFM-GTF-dominant.tre")
			exec_wqfm(input_dqrts_file, wQFM_file)

			input_dqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF.udqrts")
			wQFM_file = os.path.join(stree_output_folder, f"R{replicate_num}-wQFM-GTF-dominant-unweighted.tre")
			exec_wqfm(input_dqrts_file, wQFM_file)
			continue
		elif bucky_boot:
			input_wqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF-bucky-boot.wqrts")
			wQFM_file = os.path.join(stree_output_folder, f"R{replicate_num}-wQFM-GTF-bucky-boot.tre")
		elif boot:
			if adjusted:
				input_wqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF-boot.awqrts")
				wQFM_file = os.path.join(stree_output_folder, f"R{replicate_num}-wQFM-GTF-boot-adjusted.tre")
			else:
				input_wqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF-boot.wqrts")
				wQFM_file = os.path.join(stree_output_folder, f"R{replicate_num}-wQFM-GTF-boot.tre")
		else:
			if distance:
				if adjusted:
					input_wqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF.adwqrts")
					wQFM_file = os.path.join(stree_output_folder, f"R{replicate_num}-wQFM-GTF-distance-adjusted.tre")
				else:
					input_wqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF.dwqrts")
					wQFM_file = os.path.join(stree_output_folder, f"R{replicate_num}-wQFM-GTF-distance.tre")
			else:
				if adjusted:
					input_wqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF.awqrts")
					wQFM_file = os.path.join(stree_output_folder, f"R{replicate_num}-wQFM-GTF-adjusted.tre")
				else:
					input_wqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF.wqrts")
					wQFM_file = os.path.join(stree_output_folder, f"R{replicate_num}-wQFM-GTF.tre")
		

		exec_wqfm(input_wqrts_file, wQFM_file)
	
	with open(status_file, "a") as fout:
		fout.write("wQFM GTF done\n")

def replace_wqmc(input_wqrts_file, output_wqrts_file):
### This function converts the wqfm based quartets to wqmc compatible quartets. wqmc requires
### the taxa to be named using numbers, hence we keep a mapping to keep track. After wqmc 
### creates a tree, we replace the numbers with their corresponding taxa using the map in the
### function 'replace_back_wqmc'
	int_to_taxa = {}
	taxa_to_int = {}
	count = 1
	with open(input_wqrts_file, "r") as fin:
		for line in fin.readlines():
			line = line.split(';')[0]
			line = line.replace("(", "")
			line = line.replace(")", "")

			taxa = line.split(",")

			for taxon in taxa:
				if taxon not in taxa_to_int:
					taxa_to_int[taxon] = str(count)
					int_to_taxa[str(count)] = taxon
					count += 1

			if count > 37:
				break

	# replace now
	with open(input_wqrts_file, "r") as fin:
		with open(output_wqrts_file, "w") as fout:
			for tline in fin.readlines():
				line = tline[:]
				line = line.split(';')[0]
				line = line.replace("(", "")
				line = line.replace(")", "")

				taxa = line.split(",")
				for taxon in taxa:
					tline = tline.replace(taxon, taxa_to_int[taxon])
				tline = tline.replace("),", "|")
				tline = tline.replace("(", "")
				tline = tline.replace(")", "")
				tline = tline.replace("; ", ":")

				fout.write(tline)

	return int_to_taxa

def exec_wqmc(quartet_file, tree_file):
	cmd = "./max-cut-tree qrtt={} weights=on otre={}".format(quartet_file, tree_file)
	print(cmd)
	os.system(cmd)

def replace_back_wqmc(temporary_tree_file, output_tree_file, int_to_taxa):
### This function modifies the species tree. It restores the original name of the numbered taxa
### using the passed mapping
	# replace numbers with taxa
	with open(temporary_tree_file, "r") as fin:
		tree = fin.readline()

	for number in range(37, 0, -1):
		tree = tree.replace(str(number), int_to_taxa[str(number)])

	with open(output_tree_file, "w") as fout:
		fout.write(tree)

	os.remove(temporary_tree_file)

def run_wqmc_gtf(adjusted=False, boot=False, distance=False, dominant=False, bucky_boot=False):
	for replicate_num in range(1, REPLICATES+1):
		print(f"Replicate {replicate_num}......")
		if dominant:
			input_dqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF.dqrts")
			# translate taxa to numbers
			output_dqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF-wqmc.dqrts")
			output_tree_file = os.path.join(stree_output_folder, f"R{replicate_num}-wQMC-GTF-dominant.tre")
			# replace now
			int_to_taxa = replace_wqmc(input_dqrts_file, output_dqrts_file)
			# run wqmc
			temporary_tree_file = os.path.join(stree_output_folder, f"R{replicate_num}-wQMC-GTF-dominant.temp")
			exec_wqmc(output_dqrts_file, temporary_tree_file)
			replace_back_wqmc(temporary_tree_file, output_tree_file, int_to_taxa)

			input_dqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF.udqrts")
			# translate taxa to numbers
			output_dqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF-wqmc.udqrts")
			output_tree_file = os.path.join(stree_output_folder, f"R{replicate_num}-wQMC-GTF-dominant-unweighted.tre")
			# replace now
			int_to_taxa = replace_wqmc(input_dqrts_file, output_dqrts_file)
			# run wqmc
			temporary_tree_file = os.path.join(stree_output_folder, f"R{replicate_num}-wQMC-GTF-dominant-unweighted.temp")
			exec_wqmc(output_dqrts_file, temporary_tree_file)
			replace_back_wqmc(temporary_tree_file, output_tree_file, int_to_taxa)
			continue
		elif bucky_boot:
			input_wqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF-bucky-boot.wqrts")
			# translate taxa to numbers
			output_wqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF-bucky-boot-wqmc.wqrts")
			output_tree_file = os.path.join(stree_output_folder, f"R{replicate_num}-wQMC-GTF-bucky-boot.tre")
		elif boot:
			if adjusted:
				input_wqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF-boot.awqrts")
				# translate taxa to numbers
				output_wqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF-boot-wqmc.awqrts")
				output_tree_file = os.path.join(stree_output_folder, f"R{replicate_num}-wQMC-GTF-boot-adjusted.tre")
			else:
				input_wqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF-boot.wqrts")
				# translate taxa to numbers
				output_wqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF-boot-wqmc.wqrts")
				output_tree_file = os.path.join(stree_output_folder, f"R{replicate_num}-wQMC-GTF-boot.tre")
		else:
			if distance:
				if adjusted:
					input_wqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF.adwqrts")
					# translate taxa to numbers
					output_wqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF-wqmc.adwqrts")
					output_tree_file = os.path.join(stree_output_folder, f"R{replicate_num}-wQMC-GTF-distance-adjusted.tre")
				else:
					input_wqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF.dwqrts")
					# translate taxa to numbers
					output_wqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF-wqmc.dwqrts")
					output_tree_file = os.path.join(stree_output_folder, f"R{replicate_num}-wQMC-GTF-distance.tre")
			else:
				if adjusted:
					input_wqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF.awqrts")
					# translate taxa to numbers
					output_wqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF-wqmc.awqrts")
					output_tree_file = os.path.join(stree_output_folder, f"R{replicate_num}-wQMC-GTF-adjusted.tre")
				else:
					input_wqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF.wqrts")
					# translate taxa to numbers
					output_wqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF-wqmc.wqrts")
					output_tree_file = os.path.join(stree_output_folder, f"R{replicate_num}-wQMC-GTF.tre")

		# generate wqmc compatible quartet files
		int_to_taxa = replace_wqmc(input_wqrts_file, output_wqrts_file)

		# run wqmc		
		temporary_tree_file = os.path.join(stree_output_folder, f"R{replicate_num}-wQMC-GTF.temp")
		exec_wqmc(output_wqrts_file, temporary_tree_file)

		# replace numbers with taxa
		replace_back_wqmc(temporary_tree_file, output_tree_file, int_to_taxa)

	with open(status_file, "a") as fout:
		fout.write("wQMC GTF done\n")

def run_astral(boot=False):
	for replicate_num in range(1, REPLICATES+1):
		if boot:
			combined_gene_tree_file = os.path.join(input_folder, f'R{replicate_num}/all_boot_gt.tre')
			astral_tree_file = os.path.join(stree_output_folder, f"R{replicate_num}-astral-boot.tre")
		else:
			combined_gene_tree_file = os.path.join(input_folder, f'R{replicate_num}/all_gt.tre')
			astral_tree_file = os.path.join(stree_output_folder, f"R{replicate_num}-astral-regular.tre")

		cmd = f"java -jar {astral_path} -i {combined_gene_tree_file} -o {astral_tree_file} -t 0"
		print(cmd)
		os.system(cmd)

def compute_fp_fn_rates(additional=None):
	# suffixes = ["SVD", "wQFM-RAxML", "wQFM-SVD-exponential", "wQFM-SVD-reciprocal", "wQFM-GTF",
	# 			"wQMC-SVD-exponential", "wQMC-SVD-reciprocal", "wQMC-GTF"]

	suffixes = ["SVD", "wQFM-RAxML", "wQFM-SVD-exponential", "wQFM-SVD-reciprocal", "wQFM-GTF", "wQFM-GTF-boot", 
	     		"wQFM-GTF-adjusted", "wQFM-GTF-boot-adjusted", "wQFM-GTF-distance",
				"wQMC-SVD-exponential", "wQMC-SVD-reciprocal", "wQMC-GTF", "wQMC-GTF-boot", 
				"wQMC-GTF-adjusted", "wQMC-GTF-boot-adjusted", "wQMC-GTF-distance",
				"wQFM-GTF-dominant", "wQFM-GTF-dominant-unweighted", "wQMC-GTF-dominant", "wQMC-GTF-dominant-unweighted",
				"wQFM-GTF-bucky-boot", "wQMC-GTF-bucky-boot",
				"astral", "astral-regular", "astral-weighted", "astral-weighted-2", "astral-weighted-MB", "astral-weighted-raxml-MB"]
	

	bucky_methods = ['bucky-conc-mrbayes',
               'bucky-conc-RAxML',
               'bucky-mrbayes',
               'bucky-RAxML',
               'wQFM-bucky-mrbayes',
               'wQFM-bucky-RAxML',
               'wQMC-bucky-mrbayes',
               'wQMC-bucky-RAxML'
               ]
	
	# suffixes.extend(["a1-"+x for x in bucky_methods])
	# if cidx <= 1:
	for alpha in ["1", "1-s500"]:
		for method in bucky_methods:
			suffixes.append(f"a{alpha}-{method}")

	if additional is not None:
		suffixes.extend(additional)

	modes = [configuration]
	# suffixes = ["wQFM-SVD-exponential", "wQFM-SVD-reciprocal", "SVD", "wQMC-SVD-exponential", "wQMC-SVD-reciprocal"]
	
	for mode in modes: # per model condition
		FP_FN_FILE = os.path.join(fp_fn_folder, "FP-FN.txt")
		# Clear the file
		with open(FP_FN_FILE, mode='w') as fout:
			fout.write("")

		model_tree_file = f"{base_folder}/37-taxon/{mode}/mammalian-model-species.tre"

		for rep in range(1, REPLICATES+1): # per replicate
			for suf in suffixes: # per technique
				stree_file = os.path.join(stree_output_folder, f"R{rep}-{suf}.tre")
				
				if not os.path.exists(stree_file):
					continue

				# tree_name_in_rf = f"{mode}-{suf}/R{rep}.tre"
				tree_name_in_rf = f"R{rep}-{suf}"

				print(f"stree_file = {stree_file}, model_tree_file = {model_tree_file}")
				with open(FP_FN_FILE, mode='a') as fout:
					fout.write("{} ".format(tree_name_in_rf))

				cmd = "java -jar phylonet_v2_4.jar rf -m {} -e {} >> {}".format(model_tree_file, stree_file, FP_FN_FILE)
				
				print(cmd)
				os.system(cmd)

def compute_rf_rates():
	modes = [configuration]
	for mode in modes:
		input_file = os.path.join(fp_fn_folder, "FP-FN.txt")
		output_file = os.path.join(rf_folder, "RF.txt")
		with open(output_file, "w") as fout:
			with open(input_file, "r") as fin:
				lines = fin.readlines()
				for line in lines:
					arr = line.split(" ")
					fn = int(arr[1])
					fp = int(arr[2])
					rf = (fn+fp)/(2 * 37 - 6)
					print(arr, fn, fp, rf)
					fout.write(arr[0] + " " + str(rf) + "\n")


def run_custom_wQFM():
	for mode in ["RAxML"]:
		for replicate_num in range(1, REPLICATES+1):
			input_wqrt_file = os.path.join(output_folder, f"R{replicate_num}-{mode}.mwqrts")
			wQFM_file = os.path.join(stree_output_folder, f"R{replicate_num}-wQFM-{mode}-modified.tre")

			## Generate wQFM
			cmd = "java -jar wQFM-v1.4.jar -i {} -o {}".format(input_wqrt_file, wQFM_file)
			print(cmd)
			os.system(cmd)

	compute_fp_fn_rates([f"wQFM-{mode}-modified"])

	compute_rf_rates()

#############################################################################################

clear_status_file()

# form_fasta_text_files()

# concat_sequences()

# run_svd_wqfm()

# run_svd_wqmc()

# run_SVD_output_tree()

# concat_gene_trees()

# generate_gtf_wqrts(boot=False)

# run_wqfm_gtf(boot=False, adjusted=False)

# run_wqmc_gtf(boot=False, adjusted=False)

# run_RaXML_wqrts_wQFM()


#################### Boot Begin ############################

# concat_bootstrap_gene_trees()
# generate_gtf_wqrts(boot=True, del_boot=True)
# run_wqfm_gtf(boot=True)
# run_wqmc_gtf(boot=True)

#################### Boot END ##############################

#################### Adjustement Begin ##################### 

# for bt in  [True, False]:

# ww = 0.8

# adjust_quartets(w1=ww, w2=1-ww, boot=False)

# run_wqfm_gtf(adjusted=True, boot=False)

# run_wqmc_gtf(adjusted=True, boot=False)

#################### Adjustement End #####################

#### Distance begin ######
# generate_distance_based_gtf_wqrts(boot=False)

# run_wqfm_gtf(boot=False, distance=True)

# run_wqmc_gtf(boot=False, distance=True)

#### Distance end ######
#################### Dominant Begin #####################
# generate_dominant_quartets()
# run_wqfm_gtf(dominant=True)
# run_wqmc_gtf(dominant=True)

#################### Dominant End #####################

#################### Bucky boot Begin #####################

# concat_bucky_boot_gene_trees()

# generate_gtf_wqrts(bucky_boot=True, del_boot=True)

# run_wqfm_gtf(bucky_boot=True)

# run_wqmc_gtf(bucky_boot=True)

#################### Bucky boot End #####################

#################### ASTRAL Begin #####################

# run_astral(boot=False)
# run_astral(boot=True)

#################### ASTRAL End #####################

compute_fp_fn_rates()

compute_rf_rates()

# run_custom_wQFM()

