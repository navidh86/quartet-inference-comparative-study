import os

####################### GLOBAL VARIABLES ########################
taxa = 10 # 10, 15, 37
configurations = {}
configurations[10] = ['lower-ILS', 'higher-ILS'] # 10
configurations[15] = ['100gene-100bp', '100gene-1000bp', '1000gene-100bp', '1000gene-1000bp'][0] # 15
configurations[37] = ['1X-200-250', '0.5X-200-500', '1X-800-1000'] # 37

configuration = configurations[taxa][0]

REPLICATES = 20
GENES = 200

####################### FOLDERS ########################
input_folder = f'/home/navidh86/Desktop/comparative_study/Alignments-WeightedQuartets-SpeciesTree-main/{taxa}-taxon/{configuration}'
output_folder = f'/home/navidh86/Desktop/comparative_study/Alignments-WeightedQuartets-SpeciesTree-main/{taxa}-taxon/Output/{configuration}'
stree_output_folder = f"/home/navidh86/Desktop/comparative_study/Alignments-WeightedQuartets-SpeciesTree-main/Estimated-Species-Trees/{taxa}-taxa/{configuration}"
fp_fn_folder = f"/home/navidh86/Desktop/comparative_study/Alignments-WeightedQuartets-SpeciesTree-main/Estimated-Species-Trees/FP_FN/{taxa}-taxa"

# create the folders if not available
os.makedirs(input_folder, exist_ok=True)
os.makedirs(output_folder, exist_ok=True)
os.makedirs(stree_output_folder, exist_ok=True)
os.makedirs(fp_fn_folder, exist_ok=True)


status_file = "status_15.txt"

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
		for replicate_num in range(1, REPLICATES+1):
			sequence_file = os.path.join(output_folder, f"R{replicate_num}-Sequences-Concatenated.nex")
			output_wqrt_file = os.path.join(output_folder, f"R{replicate_num}-{mode}.wqrts")
			wQFM_file = os.path.join(stree_output_folder, f"R{replicate_num}-wQFM-RAxML.tre")
			
			
			## Run RAxML
			cmd = "python3 raxml_runner.py {} {}".format(sequence_file, output_wqrt_file)
			print(cmd)
			os.system(cmd)

			## Generate wQFM
			cmd = "java -jar wQFM-v1.3.jar -i {} -o {}".format(output_wqrt_file, wQFM_file)
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
	

def generate_gtf_wqrts():
	for replicate_num in range(1, REPLICATES+1):
		print(f"Replicate {replicate_num}......")
		combined_gene_tree_file = os.path.join(input_folder, f'R{replicate_num}/all_gt.tre')
		output_wqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF.wqrts")

		cmd = f"python3 generate_wqrts.py {combined_gene_tree_file} {output_wqrts_file}"
		print(cmd)
		os.system(cmd)

	with open(status_file, "a") as fout:
		fout.write("generating gtf qrts done\n")

def run_wqfm_gtf():
	for replicate_num in range(1, REPLICATES+1):
		print(f"Replicate {replicate_num}......")
		input_wqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF.wqrts")
		wQFM_file = os.path.join(stree_output_folder, f"R{replicate_num}-wQFM-GTF.tre")

		cmd = f"java -jar wQFM-v1.4.jar -i {input_wqrts_file} -o {wQFM_file}"
		print(cmd)
		os.system(cmd)
	
	with open(status_file, "a") as fout:
		fout.write("wQFM GTF done\n")

def run_wqmc_gtf():
	for replicate_num in range(1, REPLICATES+1):
		print(f"Replicate {replicate_num}......")
		input_wqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF.wqrts")
		wQFM_file = os.path.join(stree_output_folder, f"R{replicate_num}-wQMC-GTF.tre")

		# translate taxa to numbers
		output_wqrts_file = os.path.join(output_folder, f"R{replicate_num}-GTF-wqmc.wqrts")
		int_to_taxa = {}
		taxa_to_int = {}
		count = 0
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

				if count == 37:
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

		# run wqmc
		output_tree_file = os.path.join(stree_output_folder, f"R{replicate_num}-wQMC-GTF.tre")
		temporary_tree_file = os.path.join(stree_output_folder, f"R{replicate_num}-wQMC-GTF.temp")
		cmd = "./max-cut-tree qrtt={} weights=on otre={}".format(output_wqrts_file, temporary_tree_file)
		print(cmd)
		os.system(cmd)

		# replace numbers with taxa
		with open(temporary_tree_file, "r") as fin:
			tree = fin.readline()

		for number in range(36, -1, -1):
			tree = tree.replace(str(number), int_to_taxa[str(number)])

		with open(output_tree_file, "w") as fout:
			fout.write(tree)

		os.remove(temporary_tree_file)

	with open(status_file, "a") as fout:
		fout.write("wQMC GTF done\n")
				


def compute_fp_fn_rates():
	# suffixes = ["SVD", "wQFM-RAxML", "wQFM-SVD-exponential", "wQFM-SVD-reciprocal", "wQFM-GTF",
	# 			"wQMC-SVD-exponential", "wQMC-SVD-reciprocal", "wQMC-GTF"]

	suffixes = ["SVD", "wQFM-SVD-exponential", "wQFM-SVD-reciprocal", "wQFM-GTF",
				"wQMC-SVD-exponential", "wQMC-SVD-reciprocal", "wQMC-GTF"]

	# modes = ["1X-800-1000"]
	# modes = ["0.5X-200-500"]
	modes = [configuration]
	# suffixes = ["wQFM-SVD-exponential", "wQFM-SVD-reciprocal", "SVD", "wQMC-SVD-exponential", "wQMC-SVD-reciprocal"]
	
	for mode in modes: # per model condition
		FP_FN_FILE = os.path.join(fp_fn_folder, f"FP-FN-37tax-{mode}.txt")
		# Clear the file
		with open(FP_FN_FILE, mode='w') as fout:
			fout.write("")

		model_tree_file = f"/home/navidh86/Desktop/comparative_study/Alignments-WeightedQuartets-SpeciesTree-main/37-taxon/{mode}/mammalian-model-species.tre"

		for rep in range(1, REPLICATES+1): # per replicate
			for suf in suffixes: # per technique
				stree_file = os.path.join(stree_output_folder, f"R{rep}-{suf}.tre")				
				tree_name_in_rf = f"{mode}-{suf}/R{rep}.tre"

				print(f"stree_file = {stree_file}, model_tree_file = {model_tree_file}")
				with open(FP_FN_FILE, mode='a') as fout:
					fout.write("{} ".format(tree_name_in_rf))

				cmd = "java -jar phylonet_v2_4.jar rf -m {} -e {} >> {}".format(model_tree_file, stree_file, FP_FN_FILE)
				
				print(cmd)
				os.system(cmd)

def compute_rf_rates():
	# modes = ["1X-800-1000"]
	# modes = ["0.5X-200-500"]
	modes = [configuration]
	for mode in modes:
		input_file = f"/home/navidh86/Desktop/comparative_study/Alignments-WeightedQuartets-SpeciesTree-main/Estimated-Species-Trees/FP_FN/37-taxa/FP-FN-37tax-{mode}.txt"
		output_file = f"/home/navidh86/Desktop/comparative_study/Alignments-WeightedQuartets-SpeciesTree-main/RF-rates/37-taxon/RF-37tax-{mode}.txt"
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


#############################################################################################

# form_fasta_text_files()

# concat_sequences()

# run_svd_wqfm()

# run_svd_wqmc()

# run_SVD_output_tree()

# # run_RaXML_wqrts_wQFM()

# concat_gene_trees()

# generate_gtf_wqrts()

# run_wqfm_gtf()

# run_wqmc_gtf()

compute_fp_fn_rates()

compute_rf_rates()




