from Bio import SeqIO
import argparse


def main(args):
	## Get inputs
	input_file = args.input_file
	output_file = args.output_file
	mode1 = args.mode1
	mode2 = args.mode2

	## Process and write to output.
	## "nexus", "fasta", "phylip"
	# SeqIO.convert("output_testing.fasta", "fasta", "output_testing.nexus", "nexus", "DNA")
	SeqIO.convert(input_file, mode1, output_file, mode2, "DNA")


if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument("-i", "--input_file", type=str, help="Input file", required=True)

	parser.add_argument("-o", "--output_file", type=str, help="Output file", required=True)

	parser.add_argument("-m1", "--mode1", type=str, help="Input's mode", required=True)

	parser.add_argument("-m2", "--mode2", type=str, help="Output's mode", required=True)


	main(parser.parse_args())
