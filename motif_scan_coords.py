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
	aa_file = "../../population_alignment/arcto_hfib_translation.fasta"
	cds_file = "../../population_alignment/arcto_hfib_CDS.fasta"

	aa_names, aa_seqs = read_fasta(aa_file)
	cds_names, cds_seqs = read_fasta(cds_file)

	seq_to_pull = argv[1] # name of sequence to pull from, format: individual_allele, i.e. 1_1
	coordinates = argv[2] # coordinates of motif (aa), format: start-stop, i.e. 300-400

	aa_start = int(coordinates.split("-")[0])
	aa_stop = int(coordinates.split("-")[1])

	cds_start = aa_start * 3
	cds_stop = aa_stop * 3

	aa_motif = pull_motif(seq_to_pull, aa_names, aa_seqs, aa_start, aa_stop)
	cds_motif = pull_motif(seq_to_pull, cds_names, cds_seqs, cds_start, cds_stop)

	print(f"AA motif:\n  Length: {len(aa_motif)}\n  {aa_motif}")
	print(f"CDS motif:\n  Length: {len(cds_motif)}\n  {cds_motif}")

	# Output a table with percent ID data to use when plotting in R
	with open("ridgeline_data.csv", "w") as file:
	# with open("ridgeline_data_aa.csv", "w") as file:
	# with open("ridgeline_data_supp_aa.csv", "w") as file:
		file.write("Population,Allele,Position,Percent\n")
		for i, cds_seq in enumerate(cds_seqs):
			name = cds_names[i]
			aa_seq = aa_seqs[i]
			pop = get_population(name)
			percents = slide_motif(name, cds_seq, cds_motif)
			# percents = slide_motif(name, aa_seq, aa_motif)
			for i, percent in enumerate(percents):
				file.write(f"{pop},{name},{i},{percent}\n")



