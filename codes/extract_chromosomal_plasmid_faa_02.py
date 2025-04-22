from Bio import SeqIO

# Input and output FASTA files
input_fasta = "reference.faa"
plasmid_output = "plasmid_proteins.faa"
non_plasmid_output = "chromosomal_proteins.faa"

# Counters
total = 0
plasmid_count = 0
non_plasmid_count = 0

# Open both output files
with open(plasmid_output, "w") as plasmid_handle, open(non_plasmid_output, "w") as non_plasmid_handle:
    for record in SeqIO.parse(input_fasta, "fasta"):
        total += 1
        if "plasmid" in record.description.lower():
            plasmid_count += 1
            SeqIO.write(record, plasmid_handle, "fasta")
        else:
            non_plasmid_count += 1
            SeqIO.write(record, non_plasmid_handle, "fasta")

# Print summary
print(f"Total sequences in input file        : {total}")
print(f"Plasmid-encoded proteins extracted   : {plasmid_count}")
print(f"Chromosomal proteins extracted       : {non_plasmid_count}")
print(f"Plasmid output written to            : {plasmid_output}")
print(f"Non-plasmid output written to        : {non_plasmid_output}")