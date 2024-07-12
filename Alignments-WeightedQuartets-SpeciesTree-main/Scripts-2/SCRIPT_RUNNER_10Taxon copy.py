import os

#input_folder = '10-taxon/higher-ILS/estimated-genetrees'
#output_folder = 'highILS-quartets' # lowILS-quartets
#stree_output_folder = os.path.join("strees-estimated", "highILS")

input_folder = '10-taxon/lower-ILS/estimated-genetrees'
output_folder = 'lowILS-quartets' # lowILS-quartets
stree_output_folder = os.path.join("strees-estimated", "lowILS")
#os.makedirs(output_folder, exist_ok=True)
#os.makedirs(stree_output_folder, exist_ok=True)

get_text_file_fastas = lambda replicate_num: "R{}-Fasta-Files.txt".format(replicate_num)
get_text_file_path_fastas = lambda replicate_num: os.path.join(output_folder, get_text_file_fastas(replicate_num))
get_sequence_file_name_fasta = lambda replicate_num: "R{}-Sequences-Concatenated.fasta".format(replicate_num)
get_sequence_file_name_nexus = lambda replicate_num: "R{}-Sequences-Concatenated.nex".format(replicate_num)

get_wqrts_file = lambda replicate_num, mode: "R{}-{}.wqrts".format(replicate_num, mode)

get_wQFM_output_file_path = lambda replicate_num, mode: os.path.join(stree_output_folder, "R{}-wQFM-{}.tre".format(replicate_num, mode))

# Num replicates: R1 -> R20
# Folders: 0001 -> 0200

def convert_to_fasta():
	for replicate_num in range(1, 20+1):
		replicate_folder = os.path.join(input_folder, "R{}".format(replicate_num))
		
		print(f"Running for {replicate_folder}")

		for folder in os.listdir(replicate_folder):

			folder_path = os.path.join(replicate_folder, folder)
			phylip_file = os.path.join(folder_path, "truegene.phy")
			fasta_file = os.path.join(folder_path, "truegene.fasta")

			cmd = "python3 convert_phylip_to_fasta_manually.py -i {} -o {}".format(phylip_file, fasta_file)

			print(cmd)
			os.system(cmd)


## To rename. Remove leading zeros. 
# https://stackoverflow.com/questions/13142347/how-to-remove-leading-and-trailing-zeros-in-a-string-python
## string.lstrip("0")
def rename_by_removing_zeros():
	for replicate_num in range(1, 20+1):
		replicate_folder = os.path.join(input_folder, "R{}".format(replicate_num))
		
		print(f"Running for {replicate_folder}")

		for folder in os.listdir(replicate_folder):
			new_folder_name = folder.lstrip("0")

			folder_path = os.path.join(replicate_folder, folder)
			new_folder_path = os.path.join(replicate_folder, new_folder_name)

			cmd = "mv {} {}".format(folder_path, new_folder_path)

			print(cmd)
			os.system(cmd)	




def form_fasta_text_files():
	for replicate_num in range(1, 21):
		replicate_folder = os.path.join(input_folder, "R{}".format(replicate_num))
		text_file_path = get_text_file_path_fastas(replicate_num)

		print(f"Running for {replicate_folder}, text_file_path = {text_file_path}")
		
	

		with open(text_file_path, mode='w') as fout:
			for folder_itr in range(1, 200+1):
				folder = "{}".format(folder_itr)
				fasta_file = os.path.join(replicate_folder, folder, "truegene.fasta")

				fout.write(fasta_file)
				fout.write("\n")

		# break


def concat_sequences_generate_SVD_run_wQFM():
	for replicate_num in range(1, 20+1):
		print(f"Running for R{replicate_num}")
		text_file_path = get_text_file_path_fastas(replicate_num)
		output_file_sequences_fasta = os.path.join(output_folder, get_sequence_file_name_fasta(replicate_num))
		output_file_sequences_nexus = os.path.join(output_folder, get_sequence_file_name_nexus(replicate_num))
		
		###### CONCATENATE SEQUENCES

		cmd = "python3 concatenate_sequences.py -i {} -o {}".format(text_file_path, output_file_sequences_fasta)  
		## -n NOT_GIVEN_FOR_NOW
		#print(cmd)
		#os.system(cmd)

		##### CONVERT TO NEXUS FORMAT
		cmd = "python3 convert_formats.py -i {} -o {} -m1 fasta -m2 nexus".format(output_file_sequences_fasta, output_file_sequences_nexus)
		#print(cmd)
		#os.system(cmd)

		###### RUN SVD-reciprocal
		mode = "SVD-reciprocal" # "SVD-reciprocal"

		output_file_wqrts = os.path.join(output_folder, get_wqrts_file(replicate_num, mode))
		cmd = "python3 generate_SVD_wqrts.py -i {} -o {} -m reciprocal".format(output_file_sequences_nexus, output_file_wqrts)

		#print(cmd)
		#os.system(cmd)

		wQFM_file = get_wQFM_output_file_path(replicate_num, mode) 		####### Run wQFM
		cmd = "java -jar wQFM-v1.3.jar -i {} -o {}".format(output_file_wqrts, wQFM_file)

		print(cmd)
		os.system(cmd)

		###### RUN SVD-exponential
		mode = "SVD-exponential" # "SVD-reciprocal"

		output_file_wqrts = os.path.join(output_folder, get_wqrts_file(replicate_num, mode))
		cmd = "python3 generate_SVD_wqrts.py -i {} -o {} -m exponential".format(output_file_sequences_nexus, output_file_wqrts)

		#print(cmd)
		#os.system(cmd)

		wQFM_file = get_wQFM_output_file_path(replicate_num, mode) 		####### Run wQFM
		cmd = "java -jar wQFM-v1.3.jar -i {} -o {}".format(output_file_wqrts, wQFM_file)


		print(cmd)
		os.system(cmd)


def run_SVD_output_tree():
	for mode in ["higher-ILS", "lower-ILS"]:
		output_folder = os.path.join('strees-estimated', "{}".format(mode))
		os.makedirs(output_folder, exist_ok=True)
		#input_folder = "ZZZZDONE-{}-quartets".format(mode)
		input_folder = "{}-quartets".format(mode)
		for rep in range(1, 21):
			sequence_file = os.path.join(input_folder, "R{}-Sequences-Concatenated.nex".format(rep))
			output_tree_file = os.path.join(output_folder, "R{}-SVD.tre".format(rep))

			cmd = "python3 generate_SVD_tree.py -i {} -o {}".format(sequence_file, output_tree_file)
			print(cmd)
			os.system(cmd)


def run_RaXML_wqrts_wQFM():
	for mode in ["higher-ILS", "lower-ILS"]:
		#wqrt_folder = sequence_folder = "ZZZZDONE-{}-quartets".format(mode)
		wqrt_folder = "{}-quartets".format(mode)
		sequence_folder = "{}-quartets".format(mode)
		stree_folder = os.path.join('strees-estimated', "{}".format(mode))

		for rep in range(1, 21):
			sequence_file = os.path.join(sequence_folder, "R{}-Sequences-Concatenated.nex".format(rep))
			output_wqrt_file = os.path.join(wqrt_folder, "R{}-RAxML.wqrts".format(rep))
			wQFM_file = os.path.join(stree_folder, "R{}-wQFM-RAxML.tre".format(rep))
			
			
			## Run RAxML
			cmd = "python3 raxml_runner.py {} {}".format(sequence_file, output_wqrt_file)
			print(cmd)
			os.system(cmd)

			## Generate wQFM
			cmd = "java -jar wQFM-v1.3.jar -i {} -o {}".format(output_wqrt_file, wQFM_file)
			print(cmd)
			os.system(cmd)
			os.remove("RAxML_info.TEST")


def compute_rf_rates():

	
	mapping = {
		"higher-ILS": "higher-ILS",
		"lower-ILS": "lower-ILS"
	}

	get_tree_name_in_rf = lambda rep, suf, mode: "{}-{}/R{}.tre".format(mode, suf, rep)
	get_est_species_tree = lambda rep, suf: "R{}-{}.tre".format(rep, suf)
	get_model_tree = lambda rep, mode: os.path.join("10-taxon", mapping.get(mode), "true-speciestrees", "R{}.true.tre".format(rep))
	#get_model_tree = lambda rep, mode: os.path.join("10-taxon", mapping.get(mode), "true-speciestrees", "R{}.true.tre.stripped".format(rep))

	# R1-SVD.tre, R1-wQFM.tre, R1-wQFM-RAxML.tre, R1-wQFM-SVD-exponential.tre, R1-wQFM-SVD-reciprocal.tre
	
	suffixes = ["SVD", "wQFM-RAxML", "wQFM-SVD-exponential", "wQFM-SVD-reciprocal", "wQFM-GTF",
				"wQMC-SVD-exponential", "wQMC-SVD-reciprocal", "wQMC-GTF"]
	modes = ["higher-ILS", "lower-ILS"]
	#suffixes = ["SVD"]
	
	for mode in modes: # per model condition
		FP_FN_FILE = "/mnt/nvme/1605005-1605006/msc_project/Alignments-WeightedQuartets-SpeciesTree-main/Estimated-Species-Trees/FP_FN/FP-FN-10tax-"+mode+".txt"
		# Clear the file
		with open(FP_FN_FILE, mode='w') as fout:
			fout.write("")
		
		stree_folder = os.path.join('Estimated-Species-Trees/10-taxa-species-trees-estimated', "{}".format(mode))
		for rep in range(1, 21): # per replicate
			for suf in suffixes: # per technique
				stree_file = os.path.join(stree_folder, get_est_species_tree(rep, suf))
				model_tree_file = get_model_tree(rep, mode)
				tree_name_in_rf = get_tree_name_in_rf(rep, suf, mode)

				# print(f"stree_file = {stree_file}, model_tree_file = {model_tree_file}")
				with open(FP_FN_FILE, mode='a') as fout:
					fout.write("{} ".format(tree_name_in_rf))
				cmd = "java -jar Scripts-2/phylonet_v2_4.jar rf -m {} -e {} >> {}".format(model_tree_file, stree_file, FP_FN_FILE)
				
				print(cmd)
				os.system(cmd)


#############################################################################################

#convert_to_fasta()

#rename_by_removing_zeros()

#form_fasta_text_files()

#concat_sequences_generate_SVD_run_wQFM()

#run_SVD_output_tree()

#run_RaXML_wqrts_wQFM()

compute_rf_rates()




