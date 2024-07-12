import os
import re
import argparse


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file_name", dest="input_file", help="Input File Name DNA Sequences in NEXUS Format", type=str, required=True)
parser.add_argument("-o", "--output_file_name", dest="output_file", help="Output File Name Weighted Quartets for wQFM", type=str, required=True)
parser.add_argument("-m", "--mode", dest="mode_wqrts", help="Mode for wqrt generation (exponential/reciprocal)", type=str, required=False, default="reciprocal")
args = parser.parse_args()

input_file = args.input_file
output_file = args.output_file
mode_wqrts = args.mode_wqrts

# input_file = "anomaly_zone.nex"
# output_file = "quartets_file"

print(f"Input File: {input_file}\nOutput File: {output_file}\nMode wqrts = {mode_wqrts}")

taxa_mapping_file = "alldata.out"
batch_command_file = "batch_svdquartets_command_generation.txt"


# Form batch command
command = '''
begin paup;
	
	exec {};
	svdq evalq=all bootstrap=no qfile={} qweights={} showscores=yes qformat=qmc replace=yes;
	savetrees file={} replace=yes;
	quit;

end;
'''.format(input_file, output_file, mode_wqrts, taxa_mapping_file)

### Write the batch command to a text file
with open(batch_command_file, mode='w') as fout:
	fout.write(command)
	fout.write("\n")

### Form the command
command = "./paup4a168_ubuntu64 {}".format(batch_command_file)

### Execute the command
os.system(command)


############## Convert to wQFM format ###############

### Get taxa naming map from taxa_mapping file
with open(taxa_mapping_file, mode='r') as fin:
	lines = [l.strip() for l in fin.readlines()]


## ----------------------------------------------- ##


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


mapping_str_int = {value:key for key,value in mapping_int_str.items()}

### Convert each quartet into taxa name
### Make Newick format suitable for wQFM

def get_quartet(line):
	# EXAMPLE: 1,2|3,4:0.7803820372
	# Replace everything with COMMAs for easier splitting
	line = line.replace("|", ",").replace(":", ",")
	return line.split(",")


def get_replaced_quartet_newick(quartet):
	taxa_replaced = [mapping_int_str[int(taxa_int)] for taxa_int in quartet[0:4]]
	qrt_newick = "(({},{}),({},{})); {}".format(taxa_replaced[0], taxa_replaced[1],\
							taxa_replaced[2], taxa_replaced[3],\
							quartet[-1]) # last is weight

	return qrt_newick


with open(output_file, mode='r') as fin:
	lines = [l.strip() for l in fin.readlines()]

with open(output_file, mode='w') as fout:
	for line in lines:
		quartet = get_quartet(line)
		quartet_newick_str = get_replaced_quartet_newick(quartet)
		# print(quartet_newick_str)
		fout.write(quartet_newick_str)
		fout.write("\n")

### Write to quartet output file

####################################################

### Remove redundant files
#os.remove(batch_command_file)
#os.remove(taxa_mapping_file)


