from sys import argv
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde


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

def pull_motif(seq_name, names, seqs, start, stop):
	"""
	Description: Pull motif at given start and stop positions from given sequence
	Inputs: seq_name (string) - header of seq to pull from 
			names (list of strings) - headers from fasta file
			seqs (list of strings) - sequences from fasta file (in the same order as names)
			start (int) - start position of motif to pull
			stop (int) - stop position of motif to pull
	Return: motif (string) - motif pulled from seq
			
	"""
	index = names.index(seq_name)
	seq = seqs[index]
	motif = seq[start:stop]
	return motif

def percent_id(seq1, seq2):
	"""
	Description: Calculate the percent identity of two sequences
	Inputs: seq1, seq2 (strings) 
	Return: percent (float)
	"""
	matches = 0
	
	for i in range(len(seq1)):
		a = seq1[i]
		b = seq2[i]

		if a == b:
			matches += 1 

	return round(matches/len(seq1)*100,2)

def slide_motif(name, seq, motif):
	"""
	Description: Calculate the percent identity of the motif with every k-mer in the seq
	Inputs: seq, motif (string)
	Return: percents (list of floats)
	"""
	percents = []
	for i in range(len(seq)-len(motif)+1):
		window = seq[i:i+len(motif)]
		percent = percent_id(motif, window)
		percents.append(str(percent))

	return percents


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
	fasta_file = argv[1]
	seq_to_pull = argv[2] # name of sequence to pull from, format: individual_allele, i.e. 1_1
	coordinates = argv[3] # coordinates of motif, format: start-stop, i.e. 300-400

	names, seqs = read_fasta(fasta_file)
	
	start = int(coordinates.split("-")[0])
	stop = int(coordinates.split("-")[1])

	motif = pull_motif(seq_to_pull, names, seqs, start, stop)

	print(f"Motif:\n  Length: {len(motif)}\n  {motif}")

	# Output a table with percent ID data to use when plotting in R
	with open("motif_scan_coords.csv", "w") as file:
		file.write("Population,Allele,Position,Percent\n")
		for i, seq in enumerate(seqs):
			name = names[i]
			pop = get_population(name)
			percents = slide_motif(name, seq, motif)
			for i, percent in enumerate(percents):
				file.write(f"{pop},{name},{i},{percent}\n")



