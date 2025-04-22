from Bio import SeqIO
import csv
import os

# Function to determine conservative/non-conservative substitution
def get_substitution_type(original, mutated):
    conservative_groups = [
        {'A', 'V', 'L', 'I', 'M'},              # Aliphatic
        {'F', 'Y', 'W'},                        # Aromatic
        {'D', 'E'},                             # Acidic
        {'K', 'R', 'H'},                        # Basic
        {'S', 'T', 'N', 'Q'},                   # Polar uncharged
        {'G', 'P', 'C'}                         # Special cases
    ]
    for group in conservative_groups:
        if original in group and mutated in group:
            return "Conservative"
    return "Non-Conservative"

# Function to map aligned positions to original (ungapped) positions
def map_alignment_to_original_positions(query_sequence):
    position_map = []
    original_pos = 0
    for aa in query_sequence:
        if aa == "-":
            position_map.append(None)
        else:
            original_pos += 1
            position_map.append(original_pos)
    return position_map

# Main processing function
def extract_mutations(msa_folder, output_root_folder):
    mutation_folder = os.path.join(output_root_folder, "mutations")
    os.makedirs(mutation_folder, exist_ok=True)

    summary_log = []

    for msa_filename in os.listdir(msa_folder):
        if not msa_filename.endswith(".faa"):
            continue

        msa_path = os.path.join(msa_folder, msa_filename)

        with open(msa_path, "r") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
            if not records:
                continue
            aligned_sequences = [str(record.seq) for record in records]
            query_sequence = aligned_sequences[0]

        original_positions = map_alignment_to_original_positions(query_sequence)
        aa_matrix = [{} for _ in range(len(query_sequence))]

        for seq in aligned_sequences:
            for index, aa in enumerate(seq):
                aa_matrix[index][aa] = aa_matrix[index].get(aa, 0) + 1

        total_sequences = len(aligned_sequences)
        total_hits = total_sequences - 1
        mutations = []

        for index, count_dict in enumerate(aa_matrix):
            query_residue = query_sequence[index]
            aligned_pos = index + 1
            original_pos = original_positions[index]

            # Deletions
            if query_residue != "-" and "-" in count_dict:
                del_count = count_dict.get("-", 0)
                percentage = round((del_count / total_hits) * 100, 2)
                mutations.append([
                    aligned_pos, original_pos, query_residue,
                    "deletion", del_count, total_hits, f"{percentage}%", "Deletion"
                ])

            # Insertions
            if query_residue == "-":
                for aa, count in count_dict.items():
                    if aa != "-":
                        percentage = round((count / total_hits) * 100, 2)
                        mutations.append([
                            aligned_pos, "N/A", "-", f"insertion: {aa}",
                            count, total_hits, f"{percentage}%", "Insertion"
                        ])

            # Substitutions
            if query_residue != "-":
                filtered_residues = {
                    aa: count for aa, count in count_dict.items()
                    if aa != "-" and aa != query_residue
                }
                for aa, count in filtered_residues.items():
                    percentage = round((count / total_hits) * 100, 2)
                    substitution_type = get_substitution_type(query_residue, aa)
                    mutations.append([
                        aligned_pos, original_pos, query_residue, aa,
                        count, total_hits, f"{percentage}%", substitution_type
                    ])

        # Write mutation CSV only if mutations found
        if mutations:
            csv_filename = msa_filename.replace(".faa", ".csv")
            output_path = os.path.join(mutation_folder, csv_filename)
            summary_log.append([msa_filename, "Mutations found", len(mutations) , total_hits])

            with open(output_path, mode="w", newline="") as file:
                writer = csv.writer(file)
                writer.writerow([
                    "Aligned Position", "Original Position", "Query Residue",
                    "Mutation Type / Mutated Residue", "Count", "Total Sequences",
                    "Percentage", "Substitution Type"
                ])
                writer.writerows(mutations)

            print(f"Processed: {msa_filename} → mutations")
        else:
            # Log non-mutated proteins
            summary_log.append([msa_filename, "No mutations found", 0 , total_hits])

    # Write final summary
    summary_path = os.path.join(output_root_folder, "mutation_summary.csv")
    with open(summary_path, mode="w", newline="") as summary_file:
        writer = csv.writer(summary_file)
        writer.writerow(["MSA File", "Status", "Number of Mutations" , "Total_hits"])
        writer.writerows(summary_log)

    print(f"\nMutation summary written to: {summary_path}")

# Run script with user input
if __name__ == "__main__":
    msa_input_folder = input("Enter path to your folder containing MSA (.faa) files: ").strip()
    output_folder = input("Enter path for output folder to store mutation results: ").strip()
    extract_mutations(msa_input_folder, output_folder)