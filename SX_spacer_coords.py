import re

def read_fasta(fasta_file):
	"""
	Description: Read in sequences from fasta file 
	Inputs: fasta_file (string) - path to fasta file 
	Return: names (list of strings) - headers from fasta file
			seqs (list of strings) - sequences from fasta file (in the same order as names)
	"""
	names = []
	seqs = []
	seq = ""
	with open(fasta_file) as file:
		for line in file:
			if line.startswith(">"):
				names.append(line.strip().strip(">"))
				if seq != "":
					seqs.append(seq)
				seq = ""
			else:
				seq += line.strip()
		seqs.append(seq)
	
	return names, seqs

def matches(seq):
	"""
	Description: Find all SXnE blocks with regular expression 
	Inputs: seq (string)
	Return: positions (list of tuples) - start and stop positions of all SXnE blocks
	"""
	positions = []
	p = re.compile("(S.){3,5}E")
	for m in p.finditer(seq):
		positions.append((m.start(),m.end()))
	return positions

def get_population(allele):
	"""
	Description: Return which population an allele is from
	"""
	pop1 = [1,2,3,4,5,6,7,8]
	pop2 = [9,10,11,12,13,14,15,16,17,18]

	individual = int(allele.split("_")[0])
	if individual in pop1:
		return "1"
	else:
		return "2"


if __name__ == "__main__":
	fasta_file = "../../population_alignment/arcto_hfib_translation.fasta"
	names, seqs = read_fasta(fasta_file)
	width = .7
	gap = 1 - width

	# Output a table with coordinates of bars to draw when plotting in R
	with open("SX_ridgeline_data.csv", "w") as file:
		file.write("Population,Allele,Position,Point\n")
		for i, seq in enumerate(seqs):
			x = .5 + (gap/2)
			name = names[i]
			pop = get_population(name)
			positions = matches(seq)

			for i in range(1,len(positions)):
				point = positions[i][0]-positions[i-1][1]
				
				file.write(f"{pop},{name},{x},{0}\n")
				file.write(f"{pop},{name},{x},{point}\n")

				file.write(f"{pop},{name},{i},{point}\n")

				x += width

				file.write(f"{pop},{name},{x},{point}\n")
				file.write(f"{pop},{name},{x},{0}\n")

				x += gap
