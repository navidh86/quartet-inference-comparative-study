from Bio import SeqIO
import argparse


def get_fasta_equivalent(line):
	(tax, sequence) = line.split(" ")
	return ">{}\n{}".format(tax, sequence)


def main(args):
	## Get inputs
	input_file = args.input_file
	output_file = args.output_file
	

	with open(output_file, mode='w') as fout:
		with open(input_file, mode='r') as fin:
			line = fin.readline() # skip first line

			while line:
				line = fin.readline()
				line = line.strip()
				if len(line) > 0:
					# 8 CTAGCATTCAAGTACTTTCTTTGTTCAAATTTCCTTTTCGCCTCGTCAACTACCGTGATTGGCAAGTACGCTCGTCGTCGTTTCATACAGTGGAGGCCTG
					fasta_seq = get_fasta_equivalent(line)	
					fout.write(fasta_seq)
					fout.write("\n")
					# print(fasta_seq)




if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument("-i", "--input_file", type=str, help="Input file", required=True)

	parser.add_argument("-o", "--output_file", type=str, help="Output file", required=True)


	main(parser.parse_args())
