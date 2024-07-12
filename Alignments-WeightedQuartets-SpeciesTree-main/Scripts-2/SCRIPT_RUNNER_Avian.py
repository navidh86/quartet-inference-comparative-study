import os
import sys
import time
from ete3 import Tree
import csv

####################### GLOBAL VARIABLES ########################
INTRONS = 2516
REPLICATES = 1
genes = []
####################### FOLDERS ########################
base_folder = '/mnt/nvme/1605005-1605006/msc_project/Alignments-WeightedQuartets-SpeciesTree-main'
avian_folder = f'{base_folder}/Biological/Avian'
input_folder_intron = f'{avian_folder}/2516_Introns/2500orthologs'
output_folder_intron = f'{avian_folder}/Output'
output_folder_all = output_folder_intron
astral_path = f"/mnt/nvme/1605005-1605006/msc_project/Alignments-WeightedQuartets-SpeciesTree-main/Scripts-2/Astral.5.7.8/Astral/astral.5.7.8.jar"
stree_output_folder = f"{base_folder}/Estimated-Species-Trees/Avian"
gene_tree_folder = f"{avian_folder}/Gene-Trees/allgenes"
fp_fn_folder = f"{base_folder}/Estimated-Species-Trees/FP_FN/Avian"
rf_folder = f"{base_folder}/Estimated-Species-Trees/RF-rates/Avian"

# create the folders if not available
os.makedirs(input_folder_intron, exist_ok=True)
os.makedirs(output_folder_intron, exist_ok=True)
os.makedirs(output_folder_all, exist_ok=True)
os.makedirs(stree_output_folder, exist_ok=True)
os.makedirs(fp_fn_folder, exist_ok=True)
os.makedirs(rf_folder, exist_ok=True)

status_file = "status_avian.txt"

def clear_status_file():
	with open(status_file, 'w'):
		pass

def set_genes():
	global genes
	for subdirs, dirs, files in os.walk(input_folder_intron):
		for dir in dirs:
			genes.append(int(dir))
	genes.sort()


def set_genes_all():
	global genes
	genes = []
	for subdirs, dirs, files in os.walk(gene_tree_folder):
		for dir in dirs:
			genes.append(dir)
	genes.sort()


def form_fasta_text_files():
	if genes == []:
		set_genes()

	replicate_folder = input_folder_intron
	text_file_path = os.path.join(output_folder_intron, f"Fasta-Files-Intron.txt")

	print(f"Running for {replicate_folder}, text_file_path = {text_file_path}")
	
	with open(text_file_path, mode='w') as fout:
		for folder in genes:
			fasta_file = os.path.join(replicate_folder, str(folder), "sate.removed.intron.original.aligned-allgap.filtered")

			fout.write(fasta_file)
			fout.write("\n")

	with open(status_file, "a") as fout:
		fout.write("fasta files formed\n")

def concat_sequences_introns():
	print("Concatenating")
	text_file_path = os.path.join(output_folder_intron, f"Fasta-Files-Intron.txt")
	output_file_sequences_fasta = os.path.join(output_folder_intron, f"Sequences-Concatenated-Intron.fasta")
	output_file_sequences_nexus = os.path.join(output_folder_intron, f"Sequences-Concatenated-Intron.nex")
	
	###### CONCATENATE SEQUENCES
	# cmd = "python3 concatenate_sequences.py -i {} -o {}".format(text_file_path, output_file_sequences_fasta)  
	## -n NOT_GIVEN_FOR_NOW
	# print(cmd)
	# os.system(cmd)

	##### CONVERT TO NEXUS FORMAT
	## Instead of using the fasta file, we directly use the phylip file provided by the authors as
	## the lengths are mismatched in the concatenated fasta file
	file_sequence_phylip = os.path.join(output_folder_intron, "aln.phy") # as found in the dataset
	# cmd = "python3 convert_formats.py -i {} -o {} -m1 fasta -m2 nexus".format(output_file_sequences_fasta, output_file_sequences_nexus)
	cmd = "python3 convert_formats.py -i {} -o {} -m1 phylip -m2 nexus".format(file_sequence_phylip, output_file_sequences_nexus)
	print(cmd)
	os.system(cmd)

	with open(status_file, "a") as fout:
		fout.write("sequences conacatenated\n")


def convert_to_nexus_all():
	print("Converting")
	output_file_sequences_nexus = os.path.join(output_folder_intron, f"Sequences-Concatenated-All.nex")
	
	file_sequence_phylip = os.path.join(output_folder_intron, "ALL.aln.reduced") # as found in the dataset
	cmd = "python3 convert_formats.py -i {} -o {} -m1 phylip -m2 nexus".format(file_sequence_phylip, output_file_sequences_nexus)
	print(cmd)
	os.system(cmd)

	with open(status_file, "a") as fout:
		fout.write("sequences conacatenated\n")

def run_svd_wqfm():
	for mode in ["SVD-reciprocal", "SVD-exponential"]:
		input_file_sequences_nexus = os.path.join(output_folder_intron, f"Sequences-Concatenated-Intron.nex")
		output_file_wqrts = os.path.join(output_folder_intron, f"{mode}-intron.wqrts")
		cmd = "python3 generate_SVD_wqrts.py -i {} -o {} -m {}".format(input_file_sequences_nexus, output_file_wqrts, mode[4:])

		print(cmd)
		os.system(cmd)

		# RUN wQFM
		wQFM_file = os.path.join(stree_output_folder, f"wQFM-{mode}-intron.tre")
		cmd = "java -jar wQFM-v1.4.jar -i {} -o {}".format(output_file_wqrts, wQFM_file)

		print(cmd)
		os.system(cmd)

	with open(status_file, "a") as fout:
		fout.write("SVD wQFM done\n")


def run_SVD_output_tree(mode="Intron"):
	start_time = time.time()
	# mode = Intron or Exon or UCE or All
	sequence_file = os.path.join(output_folder_intron, f"Sequences-Concatenated-{mode}.nex")
	# sequence_file = os.path.join(output_folder_intron, f"Sequences-Concatenated-All.nex")
	output_tree_file = os.path.join(stree_output_folder, f"SVD-{mode.lower()}.tre")

	cmd = "python3 generate_SVD_tree.py -i {} -o {}".format(sequence_file, output_tree_file)
	print(cmd)
	os.system(cmd)

	with open(status_file, "a") as fout:
		fout.write("SVD output done\n")
		fout.write(f"Time taken: {time.time()-start_time} seconds")


def run_svd_wqmc():
	for mode in ["SVD-reciprocal", "SVD-exponential"]:
		input_file_sequences_nexus = os.path.join(output_folder_intron, f"Sequences-Concatenated-Intron.nex")
		output_file_wqrts = os.path.join(output_folder_intron, f"{mode}-wqmc-intron.wqrts")
		taxa_mapping_file = os.path.join(output_folder_intron, f"{mode}-wqmc-intron.mapping")

		cmd = "python3 generate_SVD_wqrts_wqmc.py -i {} -o {} -t {} -m {}".format(input_file_sequences_nexus, output_file_wqrts, taxa_mapping_file, mode[4:])

		print(cmd)
		os.system(cmd)

		# RUN wQMC
		output_tree_file = os.path.join(stree_output_folder, f"wQMC-{mode}-intron.tre")
		temporary_tree_file = os.path.join(stree_output_folder, f"wQMC-{mode}-intron.temp")
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
	set_genes_all()

	combined_gene_tree_file = os.path.join(gene_tree_folder, f'all_gt.tre')
	with open(combined_gene_tree_file, "w") as fout:
		for gene in genes:
			if "exon" in gene:
				gene_tree_file = os.path.join(gene_tree_folder, f'{gene}/C12-RAxML_bipartitions.final')
			else:
				gene_tree_file = os.path.join(gene_tree_folder, f'{gene}/RAxML_bipartitions.final')
			with open(gene_tree_file, "r") as fin:
				tree = fin.readline()

			fout.write(tree)
	with open(status_file, "a") as fout:
		fout.write("concat gene trees done\n")

def unzip_bootstrap_gene_trees():
	set_genes_all()

	for gene in genes:
		if "exon" in gene:
			bootstrap_file = os.path.join(gene_tree_folder, f'{gene}/C12-RAxML_bootstrap.all.gz')
		else:
			bootstrap_file = os.path.join(gene_tree_folder, f'{gene}/RAxML_bootstrap.all.gz')
		
		cmd = f"gunzip {bootstrap_file}"
		print(cmd)
		os.system(cmd)

	with open(status_file, "a") as fout:
		fout.write("unzipped\n")


def concat_bootstrap_gene_trees():
	set_genes_all()
	total = len(genes) # 14446
	
	# 150 partitions, each  of size 100 or less
	for partition in range(0, 150):
	# for partition in range(9, 10):
		start = partition * 100
		end = min(start + 100, total)
		# end = start + 1 # for time taking

		print(f"Processing gene: {genes[start]}")
		
		combined_gene_tree_file = os.path.join(gene_tree_folder, f'all_boot_gt_partition_{partition}.tre')
		# combined_gene_tree_file = os.path.join(gene_tree_folder, f'all_boot_gt_partition_time.tre')
		with open(combined_gene_tree_file, "w") as fout:
			for idx in range(start, end):
				gene = genes[idx]
				if "exon" in gene:
					gene_tree_file = os.path.join(gene_tree_folder, f'{gene}/C12-RAxML_bootstrap.all')
				else:
					gene_tree_file = os.path.join(gene_tree_folder, f'{gene}/RAxML_bootstrap.all')
				with open(gene_tree_file, "r") as fin:
					trees = fin.read()

				fout.write(trees)

	with open(status_file, "a") as fout:
		fout.write("concat BS gene trees done\n")

def concat_bucky_boot_gene_trees():
	if genes == []:
		set_genes()

	combined_gene_tree_file = os.path.join(input_folder, f'all_bucky_boot.tre')
	if os.path.exists(combined_gene_tree_file):
		os.remove(combined_gene_tree_file)
	for gene in genes:
		gene_tree_file = os.path.join(input_folder, f'all_in/{gene}.in')

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

def generate_gtf_wqrts(boot=False, bucky_boot=False, del_boot=False, partition=0):
	start_time = time.time()
	if boot:
		# combined_gene_tree_file = os.path.join(gene_tree_folder, f'all_boot_gt.tre')
		combined_gene_tree_file = os.path.join(gene_tree_folder, f'all_boot_gt_partition_{partition}.tre')
		# combined_gene_tree_file = os.path.join(gene_tree_folder, f'all_boot_gt_partition_time.tre')
		# output_wqrts_file = os.path.join(output_folder_all, f"GTF-boot-all.wqrts")
		output_wqrts_file = os.path.join(output_folder_all, f"GTF-boot-all-partition-{partition}.wqrts")
		# output_wqrts_file = os.path.join(output_folder_all, f"GTF-boot-all-partition-time.wqrts")
	elif bucky_boot:
		pass
		# combined_gene_tree_file = os.path.join(input_folder, f'all_bucky_boot.tre')
		# output_wqrts_file = os.path.join(output_folder, f"GTF-bucky-boot.wqrts")
	else:
		combined_gene_tree_file = os.path.join(gene_tree_folder, f'all_gt.tre')
		output_wqrts_file = os.path.join(output_folder_all, f"GTF-all.wqrts")

	cmd = f"python3 generate_wqrts.py {combined_gene_tree_file} {output_wqrts_file}"
	print(cmd)
	os.system(cmd)

	if boot and del_boot:
		os.remove(combined_gene_tree_file)

	with open(status_file, "a") as fout:
		fout.write(f"generating avian gtf qrts done in {time.time()-start_time} seconds\n")

def generate_distance_based_gtf_wqrts(boot=False, del_boot=False):
	for replicate_num in range(1, REPLICATES+1):
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
	input_wqrts_file = os.path.join(output_folder, f"GTF.wqrts")
	output_dqrts_file = os.path.join(output_folder, f"GTF.dqrts")
	output_udqrts_file = os.path.join(output_folder, f"GTF.udqrts")

	cmd = f"java -jar generateBestWQrts.jar {input_wqrts_file} {output_dqrts_file}"
	print(cmd)
	os.system(cmd)
	cmd = f"awk '{{print $1\" \"1}}' {output_dqrts_file} > {output_udqrts_file}"
	print(cmd)
	os.system(cmd)

	with open(status_file, "a") as fout:
		fout.write("generating dominant qrts done\n")

def adjust_quartets(w1=0.9, w2=0.1, boot=False, del_boot=False):
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
	if dominant:
		input_dqrts_file = os.path.join(output_folder, f"GTF.dqrts")
		wQFM_file = os.path.join(stree_output_folder, f"wQFM-GTF-dominant.tre")
		exec_wqfm(input_dqrts_file, wQFM_file)

		input_dqrts_file = os.path.join(output_folder, f"GTF.udqrts")
		wQFM_file = os.path.join(stree_output_folder, f"wQFM-GTF-dominant-unweighted.tre")
		exec_wqfm(input_dqrts_file, wQFM_file)
		return
	elif bucky_boot:
		input_wqrts_file = os.path.join(output_folder, f"GTF-bucky-boot.wqrts")
		wQFM_file = os.path.join(stree_output_folder, f"wQFM-GTF-bucky-boot.tre")
	elif boot:
		if adjusted:
			input_wqrts_file = os.path.join(output_folder, f"GTF-boot.awqrts")
			wQFM_file = os.path.join(stree_output_folder, f"wQFM-GTF-boot-adjusted.tre")
		else:
			input_wqrts_file = os.path.join(output_folder, f"GTF-boot.wqrts")
			wQFM_file = os.path.join(stree_output_folder, f"wQFM-GTF-boot.tre")
	else:
		if distance:
			if adjusted:
				input_wqrts_file = os.path.join(output_folder, f"GTF.adwqrts")
				wQFM_file = os.path.join(stree_output_folder, f"wQFM-GTF-distance-adjusted.tre")
			else:
				input_wqrts_file = os.path.join(output_folder, f"GTF.dwqrts")
				wQFM_file = os.path.join(stree_output_folder, f"wQFM-GTF-distance.tre")
		else:
			if adjusted:
				input_wqrts_file = os.path.join(output_folder, f"GTF.awqrts")
				wQFM_file = os.path.join(stree_output_folder, f"wQFM-GTF-adjusted.tre")
			else:
				input_wqrts_file = os.path.join(output_folder, f"GTF.wqrts")
				wQFM_file = os.path.join(stree_output_folder, f"wQFM-GTF.tre")
	

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
	if dominant:
		input_dqrts_file = os.path.join(output_folder, f"GTF.dqrts")
		# translate taxa to numbers
		output_dqrts_file = os.path.join(output_folder, f"GTF-wqmc.dqrts")
		output_tree_file = os.path.join(stree_output_folder, f"wQMC-GTF-dominant.tre")
		# replace now
		int_to_taxa = replace_wqmc(input_dqrts_file, output_dqrts_file)
		# run wqmc
		temporary_tree_file = os.path.join(stree_output_folder, f"wQMC-GTF-dominant.temp")
		exec_wqmc(output_dqrts_file, temporary_tree_file)
		replace_back_wqmc(temporary_tree_file, output_tree_file, int_to_taxa)

		input_dqrts_file = os.path.join(output_folder, f"GTF.udqrts")
		# translate taxa to numbers
		output_dqrts_file = os.path.join(output_folder, f"GTF-wqmc.udqrts")
		output_tree_file = os.path.join(stree_output_folder, f"wQMC-GTF-dominant-unweighted.tre")
		# replace now
		int_to_taxa = replace_wqmc(input_dqrts_file, output_dqrts_file)
		# run wqmc
		temporary_tree_file = os.path.join(stree_output_folder, f"wQMC-GTF-dominant-unweighted.temp")
		exec_wqmc(output_dqrts_file, temporary_tree_file)
		replace_back_wqmc(temporary_tree_file, output_tree_file, int_to_taxa)
		return
	elif bucky_boot:
		input_wqrts_file = os.path.join(output_folder, f"GTF-bucky-boot.wqrts")
		# translate taxa to numbers
		output_wqrts_file = os.path.join(output_folder, f"GTF-bucky-boot-wqmc.wqrts")
		output_tree_file = os.path.join(stree_output_folder, f"wQMC-GTF-bucky-boot.tre")
	elif boot:
		if adjusted:
			input_wqrts_file = os.path.join(output_folder, f"GTF-boot.awqrts")
			# translate taxa to numbers
			output_wqrts_file = os.path.join(output_folder, f"GTF-boot-wqmc.awqrts")
			output_tree_file = os.path.join(stree_output_folder, f"wQMC-GTF-boot-adjusted.tre")
		else:
			input_wqrts_file = os.path.join(output_folder, f"GTF-boot.wqrts")
			# translate taxa to numbers
			output_wqrts_file = os.path.join(output_folder, f"GTF-boot-wqmc.wqrts")
			output_tree_file = os.path.join(stree_output_folder, f"wQMC-GTF-boot.tre")
	else:
		if distance:
			if adjusted:
				input_wqrts_file = os.path.join(output_folder, f"GTF.adwqrts")
				# translate taxa to numbers
				output_wqrts_file = os.path.join(output_folder, f"GTF-wqmc.adwqrts")
				output_tree_file = os.path.join(stree_output_folder, f"wQMC-GTF-distance-adjusted.tre")
			else:
				input_wqrts_file = os.path.join(output_folder, f"GTF.dwqrts")
				# translate taxa to numbers
				output_wqrts_file = os.path.join(output_folder, f"GTF-wqmc.dwqrts")
				output_tree_file = os.path.join(stree_output_folder, f"wQMC-GTF-distance.tre")
		else:
			if adjusted:
				input_wqrts_file = os.path.join(output_folder, f"GTF.awqrts")
				# translate taxa to numbers
				output_wqrts_file = os.path.join(output_folder, f"GTF-wqmc.awqrts")
				output_tree_file = os.path.join(stree_output_folder, f"wQMC-GTF-adjusted.tre")
			else:
				input_wqrts_file = os.path.join(output_folder, f"GTF.wqrts")
				# translate taxa to numbers
				output_wqrts_file = os.path.join(output_folder, f"GTF-wqmc.wqrts")
				output_tree_file = os.path.join(stree_output_folder, f"wQMC-GTF.tre")

	# generate wqmc compatible quartet files
	int_to_taxa = replace_wqmc(input_wqrts_file, output_wqrts_file)

	# run wqmc		
	temporary_tree_file = os.path.join(stree_output_folder, f"wQMC-GTF.temp")
	exec_wqmc(output_wqrts_file, temporary_tree_file)

	# replace numbers with taxa
	replace_back_wqmc(temporary_tree_file, output_tree_file, int_to_taxa)

	with open(status_file, "a") as fout:
		fout.write("wQMC GTF done\n")

def run_astral(boot=False):
	if boot:
		combined_gene_tree_file = os.path.join(input_folder, f'all_boot_gt.tre')
		astral_tree_file = os.path.join(stree_output_folder, f"astral-boot.tre")
	else:
		combined_gene_tree_file = os.path.join(input_folder, f'all_gt.tre')
		astral_tree_file = os.path.join(stree_output_folder, f"astral-regular.tre")

	cmd = f"java -jar {astral_path} -i {combined_gene_tree_file} -o {astral_tree_file} -t 0"
	print(cmd)
	os.system(cmd)


def capitalize_species_trees(additional=None):
	suffixes = ["SVD", "wQFM-RAxML", "wQFM-SVD-exponential", "wQFM-SVD-reciprocal", "wQFM-GTF", "wQFM-GTF-boot", 
	     		"wQFM-GTF-adjusted", "wQFM-GTF-boot-adjusted", "wQFM-GTF-distance",
				"wQMC-SVD-exponential", "wQMC-SVD-reciprocal", "wQMC-GTF", "wQMC-GTF-boot", 
				"wQMC-GTF-adjusted", "wQMC-GTF-boot-adjusted", "wQMC-GTF-distance",
				"wQFM-GTF-dominant", "wQFM-GTF-dominant-unweighted", "wQMC-GTF-dominant", "wQMC-GTF-dominant-unweighted",
				"wQFM-GTF-bucky-boot", "wQMC-GTF-bucky-boot",
				"astral", "astral-regular"]
	bucky_methods = ['bucky-conc-mrbayes',
               'bucky-conc-RAxML',
               'bucky-mrbayes',
               'bucky-RAxML',
               'wQFM-bucky-mrbayes',
               'wQFM-bucky-RAxML',
               'wQMC-bucky-mrbayes',
               'wQMC-bucky-RAxML']
	
	for alpha in ["1", "1-s500"]:
		for method in bucky_methods:
			suffixes.append(f"a{alpha}-{method}")

	if additional is not None:
		suffixes.extend(additional)

	for suf in suffixes: # per technique
		stree_file = os.path.join(stree_output_folder, f"{suf}.tre")
		
		if not os.path.exists(stree_file):
			continue

		with open(stree_file, "r") as fpr:
			tree = fpr.read().upper()
			print(suf)
			print(tree)

		with open(stree_file, "w") as fpw:
			fpw.write(tree)


def infer_branch_supports(mode="all"):
	if mode != "all":
		gene_tree_file = os.path.join(gene_tree_folder, f'all_gt_{mode}.tre')
	else:
		gene_tree_file = os.path.join(gene_tree_folder, f'all_gt.tre')

	suffixes = ["SVD", "wQFM-RAxML", "wQFM-SVD-exponential", "wQFM-SVD-reciprocal", "wQFM-GTF", "wQFM-GTF-boot", 
	     		"wQFM-GTF-adjusted", "wQFM-GTF-boot-adjusted", "wQFM-GTF-distance",
				"wQMC-SVD-exponential", "wQMC-SVD-reciprocal", "wQMC-GTF", "wQMC-GTF-boot", 
				"wQMC-GTF-adjusted", "wQMC-GTF-boot-adjusted", "wQMC-GTF-distance",
				"wQFM-GTF-dominant", "wQFM-GTF-dominant-unweighted", "wQMC-GTF-dominant", "wQMC-GTF-dominant-unweighted",
				"wQFM-GTF-bucky-boot", "wQMC-GTF-bucky-boot",
				"astral", "astral-regular"]
	
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

	
	for suf in suffixes: # per technique
		stree_file = os.path.join(stree_output_folder, f"{suf}-{mode}.tre")
		
		if not os.path.exists(stree_file):
			continue

		stree_file_support = os.path.join(stree_output_folder, f"{suf}-{mode}-support.tre")
		log_file = os.path.join(stree_output_folder, f"{suf}-{mode}-support.log")

		cmd = f"java -jar {astral_path} -q {stree_file} -i {gene_tree_file} -o {stree_file_support} 2> {log_file}"
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
				"astral", "astral-regular"]
	

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

	# modes = [configuration]
	# suffixes = ["wQFM-SVD-exponential", "wQFM-SVD-reciprocal", "SVD", "wQMC-SVD-exponential", "wQMC-SVD-reciprocal"]
	
	FP_FN_FILE = os.path.join(fp_fn_folder, "FP-FN.txt")
	# Clear the file
	with open(FP_FN_FILE, mode='w') as fout:
		fout.write("")

	model_tree_file = f"{base_folder}/Biological/Avian/avian-model-species.tre"

	for mode in ["intron", "all"]:
		for suf in suffixes: # per technique
			stree_file = os.path.join(stree_output_folder, f"{suf}-{mode}.tre")
			
			if not os.path.exists(stree_file):
				continue

			temp_tree_file = f"{stree_file}.tmp"
			with open(stree_file, "r") as fpr:
				tree = fpr.read()
				with open(temp_tree_file, "w") as fpw:
					fpw.write(tree.upper())

			# tree_name_in_rf = f"{mode}-{suf}/R{rep}.tre"
			tree_name_in_rf = f"{suf}-{mode}"

			print(f"stree_file = {stree_file}, model_tree_file = {model_tree_file}")
			with open(FP_FN_FILE, mode='a') as fout:
				fout.write("{} ".format(tree_name_in_rf))

			cmd = "java -jar phylonet_v2_4.jar rf -m {} -e {} >> {}".format(model_tree_file, temp_tree_file, FP_FN_FILE)
			
			print(cmd)
			os.system(cmd)
			os.remove(temp_tree_file)

	prev_dir = "../Estimated-Species-Trees/Avian/previous"
	for file_name in os.listdir(prev_dir):
		stree_file = os.path.join(prev_dir, file_name)
		temp_tree_file = f"{stree_file}.tmp"
		with open(stree_file, "r") as fpr:
			tree = fpr.read()
			with open(temp_tree_file, "w") as fpw:
				fpw.write(tree.upper())

		# tree_name_in_rf = f"{mode}-{suf}/R{rep}.tre"
		tree_name_in_rf = stree_file

		print(f"stree_file = {stree_file}, model_tree_file = {model_tree_file}")
		with open(FP_FN_FILE, mode='a') as fout:
			fout.write("{} ".format(tree_name_in_rf))

		cmd = "java -jar phylonet_v2_4.jar rf -m {} -e {} >> {}".format(model_tree_file, temp_tree_file, FP_FN_FILE)
		
		print(cmd)
		os.system(cmd)
		os.remove(temp_tree_file)


def compute_rf_rates():
	input_file = os.path.join(fp_fn_folder, "FP-FN.txt")
	output_file = os.path.join(rf_folder, "RF.txt")
	with open(output_file, "w") as fout:
		with open(input_file, "r") as fin:
			lines = fin.readlines()
			for line in lines:
				arr = line.split(" ")
				fn = int(arr[1])
				fp = int(arr[2])
				rf = (fn+fp)/(2 * 48 - 6)
				print(arr, fn, fp, rf)
				fout.write(arr[0] + " " + str(rf) + "\n")


def fix_names():
	mapping = {}
	map_file = "../Estimated-Species-Trees/Avian/common_names_final_nbh.csv"
	with open(map_file, "r") as fp:
		csvfile = csv.reader(fp)
		for line in csvfile:
			# print(line)
			mapping[line[0]] = line[1].strip()

	print(mapping)

	suffixes = ["SVD", "wQFM-RAxML", "wQFM-SVD-exponential", "wQFM-SVD-reciprocal", "wQFM-GTF", "wQFM-GTF-boot", 
	     		"wQFM-GTF-adjusted", "wQFM-GTF-boot-adjusted", "wQFM-GTF-distance",
				"wQMC-SVD-exponential", "wQMC-SVD-reciprocal", "wQMC-GTF", "wQMC-GTF-boot", 
				"wQMC-GTF-adjusted", "wQMC-GTF-boot-adjusted", "wQMC-GTF-distance",
				"wQFM-GTF-dominant", "wQFM-GTF-dominant-unweighted", "wQMC-GTF-dominant", "wQMC-GTF-dominant-unweighted",
				"wQFM-GTF-bucky-boot", "wQMC-GTF-bucky-boot",
				"astral", "astral-regular"]
	

	bucky_methods = ['bucky-conc-mrbayes',
               'bucky-conc-RAxML',
               'bucky-mrbayes',
               'bucky-RAxML',
               'wQFM-bucky-mrbayes',
               'wQFM-bucky-RAxML',
               'wQMC-bucky-mrbayes',
               'wQMC-bucky-RAxML'
               ]

	for alpha in ["1", "1-s500"]:
		for method in bucky_methods:
			suffixes.append(f"a{alpha}-{method}")

	
	FP_FN_FILE = os.path.join(fp_fn_folder, "FP-FN.txt")
	# Clear the file
	with open(FP_FN_FILE, mode='w') as fout:
		fout.write("")

	model_tree_file = f"{base_folder}/Biological/Avian/avian-model-species.tre"

	for mode in ["intron", "all", "all-support"]:
		for suf in suffixes: # per technique
			stree_file = os.path.join(stree_output_folder, f"{suf}-{mode}.tre")
			
			if not os.path.exists(stree_file):
				continue

			with open(stree_file, "r") as fp:
				tree_string = fp.read()
				
			tree = Tree(tree_string)
			leaves = []
			for leaf in tree:
				leaf_string = str(leaf)[3:]
				leaves.append(leaf_string)
				tree_string = tree_string.replace(leaf_string, mapping[leaf_string])

			print(tree_string)
			
			renamed_stree_file = os.path.join(stree_output_folder, f"{suf}-{mode}-mapped.tre")
			with open(renamed_stree_file, "w") as fp:
				fp.write(tree_string)


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

# clear_status_file()

# set_genes()

# form_fasta_text_files()

# concat_sequences()

# convert_to_nexus_all()

# run_svd_wqfm()

# run_svd_wqmc()

# run_SVD_output_tree(mode="All")

# concat_gene_trees()

# generate_gtf_wqrts(boot=False)

# run_wqfm_gtf(boot=False, adjusted=False)

# run_wqmc_gtf(boot=False, adjusted=False)

# run_RaXML_wqrts_wQFM()


#################### Boot Begin ############################

# unzip_bootstrap_gene_trees()
concat_bootstrap_gene_trees()
generate_gtf_wqrts(boot=True, del_boot=False, partition=int(sys.argv[1]))
# run_wqfm_gtf(boot=True)
# run_wqmc_gtf(boot=True)

#################### Boot END ##############################

#################### Adjustement Begin ##################### 

# for bt in  [True, False]:

# ww = 0.9

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

# compute_fp_fn_rates()

# compute_rf_rates()

# run_custom_wQFM()

# fix_names()

#################### Branch Support #####################
# capitalize_species_trees()
# infer_branch_supports("all")


# set_genes_all()