import sys
import argparse


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq






####################################


print_once = True

def print_seq_len_once(length, num_characters):
	global print_once
	if print_once:
		print(f"Sequence Length = {length}, Using num_characters = {num_characters}")
		if num_characters is None:
			print("None: Using all Characters")
	print_once = False


def get_concatenated_sequences(files, output_file, num_characters):
	taxa_names_sequences = {}
	
	for file in files:
		fasta_sequences = SeqIO.parse(open(file),'fasta')		

		for fasta in fasta_sequences:
			name, sequence = fasta.id, str(fasta.seq)
			sequence = sequence.strip()

			# print_seq_len_once(length=len(sequence), num_characters=num_characters)
			# print(f"File name: {file}, char len = {len(sequence)}, taxa name = {name}")

			if num_characters is not None:
				sequence = sequence[:num_characters]

			if name in taxa_names_sequences:
				taxa_names_sequences[name] = taxa_names_sequences[name] + sequence
			else:
				taxa_names_sequences[name] = sequence

	return taxa_names_sequences

	'''
	## LIMITS TO 60 chars per line.
	
	with open(output_file, mode='w') as fout:
		fout.write("")

	for name in taxa_names_sequences:
		sequence = taxa_names_sequences[name]

		simple_seq_r = SeqRecord(Seq(sequence), id=name, name=name, description=name)

		with open(output_file, mode='a') as fout:
			SeqIO.write(simple_seq_r, fout, 'fasta') # , width=32765
	'''


def write_sequences_output(d_sequences, output_file):
	with open(output_file, mode='w') as fout:
		for name in d_sequences:
			print(f"Name taxa: {name}, len sequqnces = {len(d_sequences[name])}")
			s = ">{}\n{}\n".format(name, d_sequences[name])
			fout.write(s)



def main(args):

	## Get inputs
	input_file = args.file_name_for_sequences_to_concat
	output_file = args.output_file_name_merged
	num_characters = args.num_characters

	if num_characters is not None:
		print(input_file, output_file, num_characters)

	## Read input files and store
	with open(input_file, mode='r') as fin:
		files = [line.strip() for line in fin.readlines()]

	## Read fasta files' sequences and concatenate.
	sequences_concatenated = get_concatenated_sequences(files, output_file, num_characters)

	## Write to output.
	write_sequences_output(d_sequences=sequences_concatenated, output_file=output_file)

	## Write prompt conclusion
	print(f"Finished running. Written to {output_file}")


if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument("-i", "--file_name_for_sequences_to_concat", type=str, help="Input file for sequences to concatenate (in fasta format)", required=True)

	parser.add_argument("-o", "--output_file_name_merged", type=str, help="Output file for concat sequences", required=True)

	parser.add_argument("-n", "--num_characters", type=int, help="Num Characters to keep (from head)", required=False)

	main(parser.parse_args())
