import os
import re
import argparse


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file_name", dest="input_file", help="Input File Name DNA Sequences in NEXUS Format", type=str, required=True)
parser.add_argument("-o", "--output_file_name", dest="output_file", help="Output File Name Tree SVD", type=str, required=True)
args = parser.parse_args()

input_file = args.input_file
SVD_file = args.output_file

# input_file = "anomaly_zone.nex"
# output_file = "quartets_file"

print(f"Input File: {input_file}\nOutput File: {SVD_file}")

taxa_mapping_file = "alldata.out"
batch_command_file = "batch_svdquartets_command_generation.txt"


# Form batch command
command = '''
begin paup;
	
	exec {};
	svdq evalq=all bootstrap=no showscores=yes treeInf=qfm replace=yes;
	savetrees file={} replace=yes;
	quit;

end;
'''.format(input_file, taxa_mapping_file)

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
print("--------------------------printing lines-------------")
SVD_generated_tree_raw = lines[-2]
idx = lines[-2].index('(')
SVD_generated_tree = SVD_generated_tree_raw[idx:]
print(SVD_generated_tree)
print("--------------------------printing lines end-------------")
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


def get_SVD_tree(output_tree):
	decodedTree =""
	i=0
	for k in range(len(output_tree)):
		if i == len(output_tree):
			break
		
		c = (output_tree[i])
		if c!='(' and c!=')' and c!=',' and c!=';':
			idx = i;
			for j in range(i+1,len(output_tree)):
				x = output_tree[j]
			
				if x=='(' or x==')' or x==',' or x==';':
					idx = j;
					break;
			
			key = int(output_tree[i:idx])
			decodedTree = decodedTree + mapping_int_str.get(key)
			i+= (idx-i);
			
		else:
			decodedTree = decodedTree + c
			i = i+1

	print(decodedTree)
	return decodedTree

SVD_tree_decoded = get_SVD_tree(SVD_generated_tree)
print(SVD_tree_decoded)

with open(SVD_file, mode = 'w') as fout:
	fout.write(SVD_tree_decoded)


### Write to quartet output file

####################################################

### Remove redundant files
os.remove(batch_command_file)
os.remove(taxa_mapping_file)


