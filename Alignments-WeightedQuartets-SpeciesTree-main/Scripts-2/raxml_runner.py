from Bio import SeqIO
import sys
import os

nexus_input_file = str(sys.argv[1])
phylip_file = nexus_input_file[0:nexus_input_file.rindex(".")] + ".phylip"
#quartet_file = nexus_input_file[0:nexus_input_file.rindex(".")] + "_RAxML.quartets"
quartet_file = str(sys.argv[2])


records = SeqIO.parse(nexus_input_file, "nexus")
count = SeqIO.write(records, phylip_file , "phylip")
print("Converted %i records" % count)

os.system('./raxmlHPC -m GTRGAMMA -s {} -f q -p 12345 -n TEST'.format(str(phylip_file)))
# os.system('./raxmlHPC -m GTRGAMMA -s {} -f q -p 12345 -T 4 -n TEST'.format(str(phylip_file))) 
#os.system(' raxmlHPC -m GTRGAMMA -s {} -f q -p ABCDE -n TEST' .format(str(phylip_file))) # not used

f = open("RAxML_quartets.TEST", "r")
lines = f.readlines()[2:]
tax_dict = {}
newick_format = ""
for line in lines:
	temp = line.count(' ')
	if temp==1: #tax name
		arr = line.split()
		tax_dict[arr[1]] = arr[0]
	elif temp==5:
		line = line.replace(':',' ')
		line = line.replace('|',' ')
		arr = line.split()
		newick_format = newick_format + "(("+ tax_dict[arr[0]] + "," + tax_dict[arr[1]] + "),(" + tax_dict[arr[2]] + "," + tax_dict[arr[3]] + ")); " + arr[4] +"\n"
		#print(newick_format)
print(newick_format)

f = open(quartet_file, "w")
f.write(newick_format)
f.close()

