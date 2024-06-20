from sys import argv
import random
import os

pop1 = [1,2,3,4,5,6,7,8]
pop2 = [9,10,11,12,13,14,15,16,17,18]

def read_data(fasta_file):
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
    
    table = []   
    next_ID = 0
    insertion_ID_count = {}
    ID_insertion = {}

    order = list(range(len(names)))

    for index in order:
        
        ref = alignments[index]

        inGap = False

        for i in range(len(ref)):
            amino = ref[i]

            if amino == "-" and inGap == False:
                
                inGap = True
                start = i

            if amino != "-" and inGap == True:
                
                inGap = False
                stop = i
                next_ID, insertion_ID_count, ID_insertion = pull_inserts(index, names, alignments, start, stop, table, next_ID, insertion_ID_count, ID_insertion)

    return table, insertion_ID_count, ID_insertion


def pull_inserts(index, names, alignments, start, stop, table, next_ID, insertion_ID_count, ID_insertion):
    
    for i in range(len(names)):
        
        if i == index:
            continue

        row = []

        name = names[i]
        seq = alignments[i]

        insertion = seq[start:stop]

        if insertion == "-" * len(insertion):
            continue

        full_insert = insertion
        amino_insert = insertion.replace("-","")

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

        if row not in table:
            table.append(row)
        
        else:
            ID_count[1] -= 1

    return next_ID, insertion_ID_count, ID_insertion


def output_tables(table, csv_file, insertion_ID_count, ID_insertion):
    
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
    
    alignment_file = "../../population_alignment/arcto_hfib_translation_aligned.fasta"
    alignment_file = argv[1]
    names, alignments = read_data(alignment_file)
    
    csv_file = "indels_all.csv"
    csv_file = argv[2]

    table, insertion_ID_count, ID_insertion = make_table(names, alignments)
    table.sort()

    output_tables(table, csv_file, insertion_ID_count, ID_insertion)
