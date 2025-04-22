import pandas as pd
from Bio import SeqIO

# === Step 1: Load the CSV file ===
input_csv = "chromosomal/mutation_results/mutation_summary.csv"  # <-- replace with your actual filename
df = pd.read_csv(input_csv)

print(df.columns.tolist())

# === Step 2: Extract ID from 'MSA_File' column ===
df['id'] = df['MSA File'].str.extract(r'MSA_(.+)\.faa')

# === Step 3: Define mutation status based on 'Number of Mutations' ===
df['mutation_status'] = df['Number of Mutations'].apply(lambda x: 'no mutation' if x == 0 else 'mutation')

# === Step 4: Load protein sequences from FASTA file ===
fasta_file = "chromosomal_proteins.faa"  # <-- replace with your actual FASTA file
seq_records = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

# === Step 5: Get sequence length for each ID ===
lengths = []
for pid in df['id']:
    if pid in seq_records:
        lengths.append(len(seq_records[pid].seq))
    else:
        lengths.append(None)  # Or handle missing IDs as needed

df['length'] = lengths

# === Step 6: Write the final output ===
final_df = df[['id', 'length', 'mutation_status']]
final_df.to_csv("distribution_curve_chromosomal.csv", index=False)