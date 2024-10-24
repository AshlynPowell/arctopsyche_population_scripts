import re
from matplotlib.ticker import MultipleLocator, NullLocator
import matplotlib.pyplot as plt
import numpy as np
from sys import argv

# THIS SCRIPT REQUIRES MATPLOTLIB AND NUMPY

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
	
	fasta_file = argv[1]
	allele = argv[2]
	pop = get_population(allele)

	names, seqs = read_fasta(fasta_file)
	
	index = names.index(allele)
	seq = seqs[index]

	positions = matches(seq)
	points = []

    # Calculate the distance between each SXnE block (stop - previous start)
	for i in range(1,len(positions)):
		point = positions[i][0]-positions[i-1][1]
		points.append(point)

	fig, ax = plt.subplots()
	fig.set_size_inches(12,5)
	bottom = np.zeros(len(points))
	ax.bar(list(range(1,len(points)+1)), points, width=0.5, bottom=bottom, color="#020080")
	ax.set_title(f"Population: {pop}  Allele: {allele}", loc="right", size=18)
	ax.set_xlabel("G-rich Spacer Number", size=18)
	ax.set_ylabel("Length (Amino Acids)", size=18)
	ax.tick_params(axis='both', which='major', labelsize=14)
	ax.xaxis.set_major_locator(MultipleLocator(10))
	ax.xaxis.set_minor_locator(MultipleLocator(5))
	ax.yaxis.set_major_locator(MultipleLocator(10))
	ax.yaxis.set_minor_locator(MultipleLocator(5))
	ax.set_ylim([0, 140])
	ax.margins(x=0.01)
	ax.spines[['right', 'top']].set_visible(False)
	plt.savefig("SX_spacer_" + allele + ".pdf", format='pdf',dpi=1200,bbox_inches='tight', pad_inches=0.25) 

	
