import os
import csv
from Bio import SeqIO
from Bio.Blast import NCBIXML

class BLASTProcessor:
    def __init__(self, xml_folder, fasta_file, output_folder, query_fasta):
        self.xml_folder = xml_folder  # Folder containing BLAST XML files
        self.fasta_file = fasta_file  # Combined FASTA of all genome proteins
        self.output_folder = output_folder  # Output folder for results
        self.query_fasta = query_fasta  # FASTA file with query proteins
        self.csv_output_folder = os.path.join(self.output_folder, "filtered_hits_csv")

        os.makedirs(self.output_folder, exist_ok=True)
        os.makedirs(self.csv_output_folder, exist_ok=True)

    def process_all_xml_files(self):
        xml_files = [f for f in os.listdir(self.xml_folder) if f.endswith(".xml")]
        if not xml_files:
            print("No XML files found.")
            return

        for i, xml_file in enumerate(xml_files, 1):
            print(f"\n => Processing file {i}/{len(xml_files)}: {xml_file}")
            self.process_xml_file(xml_file)

        print("\n All BLAST XML files processed. Filtered outputs saved to:", self.output_folder)

    def process_xml_file(self, xml_file):
        xml_path = os.path.join(self.xml_folder, xml_file)
        with open(xml_path, "r") as file:
            blast_record = NCBIXML.read(file)

        query_def = blast_record.query
        print(f"Query definition: {query_def}")  # Debug line to see the query_def
        query_length = blast_record.query_length
        hit_details = []

        for hit in blast_record.alignments:
            for hsp in hit.hsps:
                identity = hsp.identities / hsp.align_length
                coverage = hsp.align_length / query_length
                evalue = hsp.expect

                if identity >= 0.90 and coverage >= 0.70 and evalue < 1e-5:
                    hit_details.append((hit.hit_def, identity, coverage, evalue))
                    break  # Accept hit if one HSP satisfies the condition

        hit_def_list = set(hit[0] for hit in hit_details)
        self.extract_sequences(query_def, hit_def_list, xml_file)
        self.write_hit_defs_to_csv(hit_details, xml_file)

    def extract_sequences(self, query_def, hit_def_list, xml_file):
        def write_wrapped_sequence(handle, header, sequence, line_length=60):
            handle.write(f">{header}\n")
            for i in range(0, len(sequence), line_length):
                handle.write(sequence[i:i+line_length] + "\n")
        query_id = None
        query_sequence = None

        with open(self.query_fasta, "r") as query_file:
            for record in SeqIO.parse(query_file, "fasta"):
                if query_def and query_def in record.description:
                    query_id = record.id
                    query_sequence = str(record.seq)
                    break

        if not query_id:
            print(f"Query sequence not found for {query_def}")
            return

        output_file = os.path.join(self.output_folder, f"filtered_hits_{query_id}.faa")

        with open(output_file, "w") as outfile:
            # Write query sequence first with wrapped format
            write_wrapped_sequence(outfile, query_def, query_sequence)

            # Extract matching sequences from all_genomes.fasta
            with open(self.fasta_file, "r") as infile:
                write_seq = False
                for line in infile:
                    if line.startswith(">"):
                        header = line.strip()[1:]
                        write_seq = any(hit_def == header for hit_def in hit_def_list)
                    if write_seq:
                        outfile.write(line)

    def write_hit_defs_to_csv(self, hit_details, xml_file):
        base_name = os.path.splitext(xml_file)[0]
        csv_file = os.path.join(self.csv_output_folder, f"{base_name}_hit_defs.csv")
        with open(csv_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["Hit_Def", "Identity", "Coverage", "E-value"])
            for hit_def, identity, coverage, evalue in sorted(hit_details):
                writer.writerow([hit_def, f"{identity:.4f}", f"{coverage:.4f}", f"{evalue:.2e}"])

xml_folder = input("Enter the path to the folder containing BLAST XML files: ").strip()
fasta_file = input("Enter the path to the combined genome proteins FASTA file: ").strip()
output_folder = input("Enter the path where you want to save the output files: ").strip()
query_fasta = input("Enter the path to the query FASTA file: ").strip()

# Run the processor
processor = BLASTProcessor(xml_folder, fasta_file, output_folder, query_fasta)
processor.process_all_xml_files()

"""
EXPLANATION:

1. **Initialization**:
   - Takes paths for BLAST XMLs, combined FASTA of proteins, query FASTA, and output folder.

2. **process_all_xml_files**:
   - Iterates over all XML files and processes each one.

3. **process_xml_file**:
   - Parses each BLAST XML result file.
   - Checks each HSP per hit for:
     - Identity ≥ 90%
     - Query coverage ≥ 70%
     - E-value < 1e-5
   - Adds accepted hit info as a tuple: (hit_def, identity, coverage, evalue).
   - Then extracts sequences and writes hit info to CSV.

4. **extract_sequences**:
   - Uses query_def to find and write the query sequence.
   - Then parses the all-genomes FASTA file line by line.
   - Writes sequences where the header exactly matches one in the hit_def list.

5. **write_hit_defs_to_csv**:
   - Saves the list of accepted hit definitions and associated identity, coverage, and e-value into a CSV file.
   - The CSV filename corresponds to the XML filename used.
"""