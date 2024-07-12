#15_tax_qrtt_from_gt
#weighted quartet ke newick_to_convert_wqmc python diye run kore akta text file e adjusted format er wqrt nibo
#max_cut_tree script diye run kre estimated species tree pabo
#phylonet script diye rf pabo finally
import os
#os.system("cd ..")
#os.system("ls")
#exit()
#input_folder="10-taxon/Output/lowILS-quartets"
#output_folder="10-taxon/Output/lowILS-quartets-wqmcFormat"
#sptree_output_folder="10-taxon/Output/lowILS-wqmc-speciesTree"

#input_folder="10-taxon/Output/10-taxa-quartet_from_geneTree/lower-ILS"
#output_folder="10-taxon/Output/lowILS-quartets-wqmcFormat"
#sptree_output_folder="10-taxon/Output/lowILS-wqmc-speciesTree"

input_folder="10-taxon/Output/10-taxa-quartet_from_geneTree/higher-ILS"
output_folder="10-taxon/Output/highILS-quartets-wqmcFormat"
sptree_output_folder="10-taxon/Output/highILS-wqmc-speciesTree"

#os.makedirs(output_folder, exist_ok=True)
#os.makedirs(sptree_output_folder, exist_ok=True)
#get_wqrts_file = lambda replicate_num, mode: "R{}-{}.wqrts".format(replicate_num, mode)

#get_wQMC_output_file_path = lambda replicate_num, mode: os.path.join(sptree_output_folder, "R{}-wQMC-{}.tre".format(replicate_num, mode))
#mode_list=["SVD-exponential","SVD-reciprocal"]
#mode="SVD-reciprocal"
#mode="SVD-exponential"
mode="GTF"
def run_wQMC():
	for rep in range(1, 21):
				
				input_wqrt_file = os.path.join(input_folder, "R{}-{}.wqrts".format(rep,mode))
				output_wqrt_file=os.path.join(output_folder, "wQMC_R{}-{}.wqrts".format(rep,mode))
				wQMC_file = os.path.join(sptree_output_folder, "R{}-wQMC-{}.tre".format(rep,mode))
				
				
				## Run converter script
				cmd = "python3 Scripts-2/newick_to_wqmc_converter.py {} {}".format(input_wqrt_file, output_wqrt_file)
				print(cmd)
				os.system(cmd)

				## Generate wQMC
				cmd = "Scripts-2/max-cut-tree qrtt={} weights=on otre={}".format(output_wqrt_file,wQMC_file)
				print(cmd)
				os.system(cmd)
				#cmd = "java -jar wQFM-v1.3.jar -i {} -o {}".format(output_wqrt_file, wQFM_file)
		
def compute_rf():
	exit()


run_wQMC()
