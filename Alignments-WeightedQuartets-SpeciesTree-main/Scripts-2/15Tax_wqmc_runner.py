

#weighted quartet ke newick_to_convert_wqmc python diye run kore akta text file e adjusted format er wqrt nibo
#max_cut_tree script diye run kre estimated species tree pabo
#phylonet script diye rf pabo finally
import os

variant = '1000gene-100bp'


input_folder="15-taxon/Output/"+variant+"-quartets"

output_folder="15-taxon/Output/"+variant+"-quartets"
sptree_output_folder="Estimated-Species-Trees/15-taxa-species-trees-estimated/15_tax_"+variant



def letter_to_number(file):
	temp = "numbered.txt"
	
	letters = [chr(x) for x in range(ord('A'), ord('P'))]
	with open(temp, "w") as wp:
		with open(file, "r") as fp:
			for line in fp.readlines():
				for letter in letters:
					line = line.replace(letter, str(ord(letter)-ord('A')))
				wp.write(line)

	return temp

def number_to_letter(file):
	with open(file, "r") as fp:
		line = fp.readline()
	
	# numbers = [chr(x) for x in range(ord('A'), ord('P'))]
	with open(file, "w") as fp:
		for i in range(14, -1, -1):
			line = line.replace(str(i), chr(i+ord('A')))
		fp.write(line)

#get_wQMC_output_file_path = lambda replicate_num, mode: os.path.join(sptree_output_folder, "R{}-wQMC-{}.tre".format(replicate_num, mode))
#mode_list=["SVD-exponential","SVD-reciprocal"]
mode="SVD-reciprocal"
#mode="SVD-exponential"
#mode="GTF"
def run_wQMC():
	modes = ["SVD-exponential", "SVD-reciprocal", "GTF"]
	for mode in modes:
		for rep in range(1, 11):			
			input_wqrt_file = os.path.join(input_folder, "R{}-{}.wqrts".format(rep,mode))
			output_wqrt_file=os.path.join(output_folder, "R{}-wQMC-{}.wqrts".format(rep,mode))
			wQMC_file = os.path.join(sptree_output_folder, "R{}-wQMC-{}.tre".format(rep,mode))
			
			## Run converter script
			cmd = "python3 Scripts-2/newick_to_wqmc_converter.py {} {}".format(input_wqrt_file, output_wqrt_file)
			print(cmd)
			os.system(cmd)

			output_wqrt_file = letter_to_number(output_wqrt_file)

			# exit()

			## Generate wQMC
			cmd = "Scripts-2/max-cut-tree qrtt={} weights=on otre={}".format(output_wqrt_file,wQMC_file)
			print(cmd)
			os.system(cmd)

			os.remove(output_wqrt_file)

			number_to_letter(wQMC_file)
				#cmd = "java -jar wQFM-v1.3.jar -i {} -o {}".format(output_wqrt_file, wQFM_file)
		
def compute_rf():
	exit()


run_wQMC()
