from sys import argv
import os

def read_fasta(fasta_file, gene):
	"""
	Description: Read in sequences from fasta file 
	Inputs: fasta_file (string) - path to fasta file
			gene (string) 
	Return: names (list of strings) - headers from fasta file
			seqs (list of strings) - sequences from fasta file (in the same order as names)
	"""
	names = []
	seqs = []
	seq = ""
	with open(fasta_file) as file:
		for line in file:
			if line.startswith(">"):
				if gene == "resilin":
					names.append(line[7:-11]+line[-2])
				else:
					names.append(line[1:-1])
				if seq != "":
					seqs.append(seq)
				seq = ""
			else:
				seq += line.strip()
		seqs.append(seq)

	return names, seqs

def remove_common_gaps(seq1, seq2):
	"""
	Description: Remove gaps shared between two sequences (resulting in a pairwise alignment)
	Inputs: seq1 (string) 
			seq2 (string)
	Return: new1 (string) - seq1 with shared gaps removed
			new2 (string) - seq2 with shared gaps removed
	"""
	new1 = ""
	new2 = ""
	for i in range(len(seq1)):
		amino1 = seq1[i]
		amino2 = seq2[i]

		if not (amino1 == "-" and amino2 == "-"):
			new1 += amino1
			new2 += amino2

	return new1, new2

def get_indels(seq1, seq2):
	"""
	Description: Pull indels from pairwise alignment
	Inputs: seq1 (string)
			seq2 (string)
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

		if inGap == True and ((allele == "1" and amino2 != "-") or (allele == "2" and amino1 != "-")):
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

	# Input alignment file and gene, read in alignment data
	alignment_file = argv[1]
	gene = argv[2]
	names, alignments = read_fasta(alignment_file, gene)

	order = ['1_1', '1_2', '2_1', '2_2', '3_1', '3_2', '4_1', '4_2', '5_1', '5_2', '6_1', '6_2', '7_1', 
          	 '7_2', '8_1', '8_2', '9_1', '9_2', '10_1', '10_2', '11_1', '11_2', '12_1', '12_2', '13_1', 
          	 '13_2', '14_1', '14_2', '15_1', '15_2', '16_1', '16_2', '17_1', '17_2', '18_1', '18_2']

	if gene == "resilin":
		order_old = order
		order = []
		for i in ["A","B"]:
			for j in order_old:
				order.append(j+i)
	else:
		order_old = order
		order = []
		for i in order_old:
			if i in names:
				order.append(i)

	# Loop over each pair of sequences in population alignment
	for i in range(len(order)-1):

		name1 = order[i]
		name2 = order[i+1]

		seq1 = alignments[names.index(name1)]
		seq2 = alignments[names.index(name2)]

		# Path to csv file for coords of pairwise alignment between two given sequences
		csv_file = gene + "/" + name1 + "vs" + name2 + "coords.csv"

		# Remove common gaps and get positions of indels
		new1, new2 = remove_common_gaps(seq1, seq2)
		table = fix_positions(get_indels(new1, new2))
		
		# Calculate output polygon coords for plotting in R 
		coords = get_poly_coords(table, new1, new2)
		output_coords(coords, csv_file)


