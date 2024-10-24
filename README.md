# Arctopsyche Population Scripts
This repository includes the scripts used to analyze silk genes in populations of Arctopsyche grandis.

## Figure 1: Alleles of single individuals from multiple species of caddisfly

`pairwise_alignment_coords.py` performs a pairwise alignment of two alleles with MUSCLE then calculates the coordinates of each indel in the alignment and outputs the coordinates to a csv file. 

`pairwise_alignment_fig.R` uses the csv file outputted from `pairwise_alignment_coords.py` to create a plot of the pairwise alignment that highlights allele length and indel size and location.

## Figure 3: Indel histogram 

`pull_indels.py` pulls all the indels from a multiple sequence alignment and outputs the data, including indel lengths and counts, to a csv file.

`plot_indel_lengths.py` creates a histogram of the indel lengths from the csv file outputted by `pull_indels.py`. Given that most of the indels are of smaller size, the plot only includes the counts for indels up to 100 amino acids.

## Figure 4: Motif amino acid sequence scan 

`motif_scan_coords.py` utilizes a sliding window approach to calculate the percent identity of a given motif with every k-mer in each allele. The values are outputted to a csv file to be used in plotting. 

`motif_scan_fig.R` creates a ridgeline plot of the percent identities calculated with `motif_scan_coords.py`.

## Figure 5: (SX)nE spacer patterns

`SX_spacer_coords.py` calculates the distance between serine blocks in each allele and outputs them to a csv file. 

`SX_spacer_fig.R` plots all of the distances calculated with `SX_spacer_coords.py` in a ridgeline-style figure.

## Supplemental figures: motif scans and (SX)nE spacer patterns from individual alleles 

`motif_scan_single.py` and `SX_spacer_single.py` create the motif scan and SX spacer plots on individual alleles. 


