import os
import csv
from Bio import SeqIO

# === Take input from user ===
summary_csv = input("Enter path to mutation summary CSV file: ").strip()
hit_def_folder = input("Enter path to folder containing hit definition CSVs: ").strip()
combined_fasta_file = input("Enter path to reference FASTA (.faa) file: ").strip()
output_folder = input("Enter output folder path: ").strip()

# === Auto-generated output filenames ===
output_mutation_csv = os.path.join(output_folder, "mutation_summary_with_hit_defs.csv")
output_non_mutation_csv = os.path.join(output_folder, "non_mutation_summary_with_hit_defs.csv")

# === Step 1: Read FASTA headers into a dictionary: {protein_id: description} ===
reference_defs = {}
with open(combined_fasta_file) as fasta:
    for line in fasta:
        if line.startswith(">"):
            line = line.strip()
            parts = line[1:].split(maxsplit=1)
            if len(parts) == 2:
                protein_id, description = parts
            else:
                protein_id = parts[0]
                description = "Unknown"
            reference_defs[protein_id] = description

# === Step 2: Load full sequences from FASTA for length lookup ===
seq_records = SeqIO.to_dict(SeqIO.parse(combined_fasta_file, "fasta"))

# === Step 3: Containers for outputs ===
mutation_rows = []
non_mutation_rows = []

# === Step 4: Process mutation summary ===
with open(summary_csv, newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        msa_filename = row["MSA File"]
        status = row["Status"]
        num_mutations = row["Number of Mutations"]

        # Extract protein_id
        protein_id = msa_filename.replace("MSA_", "").replace(".faa", "")

        # Get Hit Def
        hit_def_path = os.path.join(hit_def_folder, f"{protein_id}_hit_defs.csv")
        if os.path.exists(hit_def_path):
            with open(hit_def_path, newline="") as hitfile:
                hit_reader = csv.reader(hitfile)
                next(hit_reader, None)  # skip header
                try:
                    hit_def = next(hit_reader)[0]
                except StopIteration:
                    hit_def = "Unknown"
        else:
            hit_def = "Unknown"

        # Get Reference Def
        reference_def = reference_defs.get(protein_id, "Unknown")

        # Get sequence length
        length = len(seq_records[protein_id].seq) if protein_id in seq_records else "Unknown"

        # Prepare output row
        record = [protein_id, reference_def, hit_def,
                  status.replace("Mutations found", "Mutated").replace("No mutations", "Not Mutated"),
                  num_mutations, length]

        if status == "Mutations found":
            mutation_rows.append(record)
        else:
            non_mutation_rows.append(record)

# === Step 5: Write outputs with length included ===
headers = ["Protein_Id", "Reference_Def", "Hit Def", "Status", "Number of Mutations", "Length"]

with open(output_mutation_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(headers)
    writer.writerows(mutation_rows)

with open(output_non_mutation_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(headers)
    writer.writerows(non_mutation_rows)

print("  →", output_mutation_csv)
print("  →", output_non_mutation_csv)