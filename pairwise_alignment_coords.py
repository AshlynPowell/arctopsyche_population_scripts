from sys import argv
import os

# THIS SCRIPT REQUIRES MUSCLE TO BE INSTALLED

def pairwise_alignment(fasta_file, alignment_file):
	"""
	Description: Perform a pairwise alignment using MUSCLE 
	Inputs: fasta_file (string) - path to fasta file (used as input file for MUSCLE)
			alignment_file (string) - path to alignment fasta file (used as output file for MUSCLE)
	"""
	os.system(f"muscle -align {fasta_file} -output {alignment_file}")

def read_fasta(fasta_file):
	"""
	Description: Read in sequences from fasta file 
	Inputs: fasta_file (string) - path to fasta file 
	Return: seqs (list of strings) - sequences from fasta file
	"""
	seqs = []
	seq = ""
	with open(fasta_file) as file:
		for line in file:
			if line.startswith(">"):
				if seq != "":
					seqs.append(seq)
				seq = ""
			else:
				seq += line.strip()
		seqs.append(seq)
	
	return seqs

def get_indels(seq1, seq2):
	"""
	Description: Pull indels from pairwise alignment
	Inputs: seq1 (string) - first sequence in the alignment
			seq2 (string) - second sequence in the alignment
	Return: table (list of lists of strings) - table of indel data with columns: allele (with insert), start position, indel length
	"""
	inGap = False
	table = []

	# Loop over each position in the alignment
	for i in range(len(seq1)):
		
		# Pull amino acid for given position in each seq
		amino1 = seq1[i]
		amino2 = seq2[i]

		if (amino1 == "-" or amino2 == "-") and inGap == False:
			# This would be the start of a gap 

			inGap = True
			row = []

			# Add which allele has the insertion 
			if amino1 == "-":
				allele = "2"
				row.append(allele)
			else:
				allele = "1"
				row.append(allele)

			# Add the start position
			start = i
			row.append(str(i))


		if amino1 != "-" and amino2 != "-" and inGap == True:
			# This would be the end of a gap

			inGap = False

			# Pull the insertion from the sequence
			if allele == "1":
				insertion = seq1[start:i]
			else:
				insertion = seq2[start:i]

			# Add the length of the insertion
			row.append(str(len(insertion)))
			table.append(row)

	return table

def fix_positions(table):
	"""
	Description: Adjust the values in the table so that indel start positions are in relation to the unaligned sequences
	Inputs: table (list of lists of strings) - table of indel data with columns: allele (with insert), start position, indel length
	Return: table (list of lists of strings) - table of indel data with adjusted start positions (same columns)
	"""

	offset1 = 0
	offset2 = 0

	for i, row in enumerate(table):
		allele = row[0]
		position = int(row[1])
		length = int(row[2])

		if allele == "1":
			row[1] = position - offset1
			table[i] = row
			offset1 += length

		else:
			row[1] = position - offset2
			table[i] = row
			offset2 += length

	return table

def get_poly_coords(table, seq1, seq2):
	"""
	Description: Calculate the coordinates of polygons based on indel positions (the grey boxes in the end figure)
	Inputs: table (list of lists of strings) - table of indel data with columns: allele (with insert), start position, indel length
			seq1 (string) - first sequence in the alignment
			seq2 (string) - second sequence in the alignment
	Return: coords (list of lists of strings) - table of coordinates of polygons with columns: top_left, top_right, bottom_left, bottom_right
	"""
	if table == []:
	  return [["0",str(len(seq1)),"0",str(len(seq2))]]
	
	coords = []
	
	prev_top = 0
	prev_bottom = 0
	prev_length = 0

	offset = 0

	# Loop over each row (i.e. each indel) in the table and add coodinates of each box
	for row in table:
		allele = row[0]
		position = int(row[1])
		length = int(row[2])

		if allele == "1":
			top = position - offset
			bottom = position

		else:
			top = position
			bottom = position + offset

		new_row = [str(prev_top), str(top), str(prev_bottom), str(bottom)]
		coords.append(new_row)

		if allele == "1":
			offset -= length
			top += length

		else:
			offset += length
			bottom += length

		prev_top = top
		prev_bottom = bottom
		prev_length = length

	last_row = [str(top), str(len(seq1.replace("-",""))), str(bottom), str(len(seq2.replace("-","")))]
	coords.append(last_row)

	return coords


def output_coords(coords, csv_file):
	"""
	Description: Output polygon coordinates to a csv file
	Inputs: coords (list of lists of strings) - table of coordinates of polygons with columns: top_left, top_right, bottom_left, bottom_right
			csv_file (string) - path to csv file to output to
	"""
	with open(csv_file, "w") as file:
		header = ["X1", "X2", "X3", "X4"]
		line = ",".join(header) + "\n"
		file.write(line)
		
		for row in coords:
			line = ",".join(row) + "\n"
			file.write(line)


if __name__ == "__main__":
	
	# Input species code and create strings with paths to input/output files using the code
	code = argv[1] 
	fasta_file = "../alignments/" + code + "_both_alleles.fasta"
	alignment_file = "../alignments/" + code + "_both_alleles_aligned.fasta"
	csv_file = code + "_both_alleles.csv"
	
	# Create and read in pairwise alignment 
	pairwise_alignment(fasta_file, alignment_file)
	alignments = read_fasta(alignment_file)

	# Get positions of indels
	table = get_indels(alignments[0], alignments[1])
	table = fix_positions(table)
	
	# Calculate output polygon coords for plotting in R 
	coords = get_poly_coords(table, alignments[0], alignments[1])
	output_coords(coords, csv_file)
	
