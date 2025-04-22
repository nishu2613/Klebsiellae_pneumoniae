import os
from Bio import SeqIO
import subprocess
from concurrent.futures import ThreadPoolExecutor

# ---- User Inputs ----
protein_type = input("Which proteins do you want to BLAST? (plasmid/chromosomal): ").strip().lower()

# Validate input
while protein_type not in ["plasmid", "chromosomal"]:
    protein_type = input("Please enter 'plasmid' or 'chromosomal': ").strip().lower()

query_fasta = input(f"Enter your {protein_type} protein FASTA file: ").strip()
blast_db = input("Enter BLAST database path: ").strip()
output_dir = input(f"Enter output directory for {protein_type} XML results: ").strip()
threads = int(input("Enter number of threads to use (e.g. 12): ").strip())

# ---- Setup Output Folder ----
os.makedirs(output_dir, exist_ok=True)
temp_fasta_dir = os.path.join(output_dir, f"temp_fastas_{protein_type}")
os.makedirs(temp_fasta_dir, exist_ok=True)

# ---- Prepare Individual FASTA Files ----
print(f"\nSplitting {protein_type} protein sequences into individual FASTA files...")
fasta_paths = []

for record in SeqIO.parse(query_fasta, "fasta"):
    fasta_path = os.path.join(temp_fasta_dir, f"{record.id}.fasta")
    SeqIO.write(record, fasta_path, "fasta")
    fasta_paths.append(fasta_path)

print(f"Total sequences to BLAST: {len(fasta_paths)}")

# ---- BLAST Function ----
def run_blast(fasta_file):
    base_name = os.path.splitext(os.path.basename(fasta_file))[0]
    xml_output = os.path.join(output_dir, f"{base_name}.xml")

    cmd = [
        "blastp",
        "-query", fasta_file,
        "-db", blast_db,
        "-out", xml_output,
        "-outfmt", "5"
    ]
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

# ---- Run BLAST in Parallel ----
print(f"\nRunning BLASTP on {len(fasta_paths)} {protein_type} proteins using {threads} threads...\n")

with ThreadPoolExecutor(max_workers=threads) as executor:
    executor.map(run_blast, fasta_paths)

print(f"\nAll BLASTP results for {protein_type} proteins saved as XML in:", output_dir)