from sys import argv
import os
 
def read_data(fasta_file):
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

def pairwise_alignment(fasta_file, alignment_file):
	
	# os.system(f"mafft --auto {fasta_file} > {alignment_file}")

	os.system(f"muscle -align {fasta_file} -output {alignment_file}")

def get_indels(seq1, seq2):
	
	inGap = False
	table = []
	
	next_ID = 0
	insertion_ID_count = {}
	ID_insertion = {}

	for i in range(len(seq1)):
		amino1 = seq1[i]
		amino2 = seq2[i]

		if (amino1 == "-" or amino2 == "-") and inGap == False:
			
			inGap = True
			row = []

			if amino1 == "-":
				allele = "2"
				row.append(allele)
			else:
				allele = "1"
				row.append(allele)

			start = i
			row.append(str(i))

		if amino1 != "-" and amino2 != "-" and inGap == True:
			inGap = False

			if allele == "1":
				insertion = seq1[start:i]
			else:
				insertion = seq2[start:i]

			row.append(str(len(insertion)))
			table.append(row)

	return table

def fix_positions(table):
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


	with open(csv_file, "w") as file:
		header = ["X1", "X2", "X3", "X4"]
		line = ",".join(header) + "\n"
		file.write(line)
		
		for row in coords:
			line = ",".join(row) + "\n"
			file.write(line)


if __name__ == "__main__":
	
	code = argv[1] 
	fasta_file = "../alignments/" + code + "_both_alleles.fasta"
	alignment_file = "../alignments/" + code + "_both_alleles_aligned.fasta"
	csv_file = code + "_both_alleles.csv"
	
	pairwise_alignment(fasta_file, alignment_file)

	alignments = read_data(alignment_file)

	print(get_indels(alignments[0], alignments[1]))

	table = fix_positions(get_indels(alignments[0], alignments[1]))
	coords = get_poly_coords(table, alignments[0], alignments[1])
	output_coords(coords, csv_file)
	





