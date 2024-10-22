from sys import argv
import random
import os

pop1 = [1,2,3,4,5,6,7,8]
pop2 = [9,10,11,12,13,14,15,16,17,18]

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

def make_table(names, alignments):
    """
    Description: Make a table of indel data
    Inputs: names (list of strings) - headers from fasta file
            alignments (list of strings) - sequences from the alignment file (in the same order as names)
    Return: table (list of lists of strings) - table with indel positions and lengths with 
                columns: population,individual,allele,position,ID,full length,amino acid length
            insertion_ID_count (dictionary, key=insertion, value=[ID,count]) - keeps track of the ID and count for each insertion
            ID_insertion (dictionary, key=ID, value=insertion) - for quick look up of insertion by ID
    """
    
    table = []   
    next_ID = 0
    insertion_ID_count = {}
    ID_insertion = {}
    checked = []

    for index in range(len(names)):
        # Set reference sequence
        ref = alignments[index]

        inGap = False

        # Loop over reference sequence to identify gaps and pull indels
        for i in range(len(ref)):
            amino = ref[i]

            if amino == "-" and inGap == False:
                # This is the start of a gap
                inGap = True
                start = i

            if amino != "-" and inGap == True:
                # This is the end of a gap
                inGap = False
                stop = i
                # This skips over positions that have previously been checked 
                # (i.e. this gap was already found in previous reference sequence)
                if (start,stop) in checked:
                    continue
                
                checked.append((start,stop))

                # Pull indels from all other sequences 
                next_ID, insertion_ID_count, ID_insertion = pull_inserts(index, names, alignments, start, stop, table, next_ID, insertion_ID_count, ID_insertion)

    return table, insertion_ID_count, ID_insertion


def pull_inserts(index, names, alignments, start, stop, table, next_ID, insertion_ID_count, ID_insertion):
    """
    Description: Pull indels at given start and stop positions and add them to the table
    Inputs: index (integer) - index of the reference sequence (in names and alignments lists)
            names (list of strings) - headers from fasta file
            alignments (list of strings) - sequences from the alignment file (in the same order as names)
            start (integer) - start position of gap
            stop (integer) - stop position of gap 
            table (list of lists of strings) - table with indel positions and lengths with 
                columns: population,individual,allele,position,ID,full length,amino acid length
            next_ID (integer) - next available ID number (for new indels found)
            insertion_ID_count (dictionary, key=insertion, value=[ID,count]) - keeps track of the ID and count for each insertion
            ID_insertion (dictionary, key=ID, value=insertion) - for quick look up of insertion by ID
    Return: table - same as above, just updated with new indels 
            insertion_ID_count - same as above, just updated with new indels 
            ID_insertion - same as above, just updated with new indels 
    """    
    found = []
    
    # Loop over all sequence in the alignment, skipping the reference sequence
    for i in range(len(names)):
        if i == index:
            continue

        row = []

        name = names[i]
        seq = alignments[i]

        # Pull indel at given start and stop, skipping when indel is also a gap or previously found 
        insertion = seq[start:stop]

        if insertion == "-" * len(insertion) or insertion in found:
            continue

        # Keep a list of indels found at this location so we only count unique ones 
        found.append(insertion)

        full_insert = insertion
        amino_insert = insertion.replace("-","")

        # Add insert to dictionaries if new or update count if already existing 
        if amino_insert in insertion_ID_count:
            ID_count = insertion_ID_count[amino_insert]
            ID = ID_count[0]
            ID_count[1] += 1
            insertion_ID_count[amino_insert] = ID_count
        else:
            ID = next_ID
            insertion_ID_count[amino_insert] = [next_ID, 1]
            ID_insertion[next_ID] = amino_insert
            next_ID += 1

        # Add data to table 
        individual = int(name.split("_")[0])
        if individual in pop1:
            row.append("1")
        elif individual in pop2:
            row.append("2")
        else:
            print("ISSUE")

        row += name.split("_")
        row.append(start)
        row.append(ID)
        row.append(len(full_insert))
        row.append(len(amino_insert))
        
        row = [int(x) for x in row]

        table.append(row)

    return next_ID, insertion_ID_count, ID_insertion


def output_tables(table, csv_file, insertion_ID_count, ID_insertion):
    """
    Description: Output tables to csv file to be used in plotting 
    Inputs: table (list of lists of strings) - table with indel positions and lengths with 
                columns: population,individual,allele,position,ID,full length,amino acid length
            csv_file (string) - path to output csv file 
            insertion_ID_count (dictionary, key=insertion, value=[ID,count]) - keeps track of the ID and count for each insertion
            ID_insertion (dictionary, key=ID, value=insertion) - for quick look up of insertion by ID
    """
    with open(csv_file, "w") as file:
        header = ["Population", "Individual", "Allele", "Position", "ID", "Full Length", "Amino Acid Length"]
        line = ",".join(header) + "\n"
        file.write(line)
        
        for row in table:
            row = [str(x) for x in row]
            line = ",".join(row) + "\n"
            file.write(line)

        file.write("\n\n\n")

        header = ["ID", "Count", "Insertion"]

        line = ",".join(header) + "\n"
        file.write(line)

        for ID, insertion in ID_insertion.items():
            ID_count = insertion_ID_count[insertion]
            count = ID_count[1]
            line = f"{ID},{count},{insertion}\n"

            file.write(line)

if __name__ == "__main__":
    
    alignment_file = argv[1]
    csv_file = argv[2]
    names, alignments = read_fasta(alignment_file)

    table, insertion_ID_count, ID_insertion = make_table(names, alignments)
    table.sort()

    output_tables(table, csv_file, insertion_ID_count, ID_insertion)
