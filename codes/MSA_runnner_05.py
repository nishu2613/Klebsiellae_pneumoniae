import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from threading import Lock

# Global locks for writing safely from parallel processes
log_lock = Lock()
skipped_lock = Lock()

def count_fasta_sequences(file_path):
    with open(file_path) as f:
        return sum(1 for line in f if line.startswith(">"))


def run_single_mafft_job(input_folder, output_folder, filename, log_file, skipped_file):
    input_path = os.path.join(input_folder, filename)
    base_name = filename.replace("filtered_hits_", "")
    output_filename = f"MSA_{base_name}"
    output_path = os.path.join(output_folder, output_filename)

    # Skip files with fewer than 2 sequences
    if count_fasta_sequences(input_path) < 2:
        with skipped_lock:
            with open(skipped_file, "a") as skipped:
                skipped.write(f"{filename}\n")
        with log_lock:
            with open(log_file, "a") as log:
                log.write(f"[{filename}] Skipped: Less than 2 sequences.\n")
        return f"{filename} skipped"

    # Run MAFFT
    try:
        result = subprocess.run(
            ["mafft", "--auto", input_path],
            stdout=open(output_path, "w"),
            stderr=subprocess.PIPE,
            text=True
        )
        if result.stderr:
            with log_lock:
                with open(log_file, "a") as log:
                    log.write(f"\n[{filename}]\n{result.stderr}\n")
        return f"{filename} aligned"
    except Exception as e:
        with log_lock:
            with open(log_file, "a") as log:
                log.write(f"\n[{filename}] Failed to run MAFFT: {str(e)}\n")
        return f"{filename} failed"


def run_mafft_on_folder_parallel(input_folder, output_folder, log_file, skipped_file, max_workers=4):
    os.makedirs(output_folder, exist_ok=True)

    faa_files = [f for f in os.listdir(input_folder) if f.endswith(".faa") and f.startswith("filtered_hits_")]
    total_files = len(faa_files)

    # Clear logs
    open(log_file, "w").close()
    open(skipped_file, "w").close()

    print(f"Total files to process: {total_files}")

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(run_single_mafft_job, input_folder, output_folder, f, log_file, skipped_file)
            for f in faa_files
        ]
        for i, future in enumerate(as_completed(futures), 1):
            try:
                result = future.result()
                print(f"[{i}/{total_files}] {result}")
            except Exception as e:
                print(f"[{i}/{total_files}] Error: {str(e)}")


# Ask user for input
if __name__ == "__main__":
    input_folder = input("Enter the input folder name: ").strip()
    output_folder = input("Enter the output folder name: ").strip()
    log_file = input("Enter the name for the log file (e.g., mafft_errors.log): ").strip()
    skipped_file = input("Enter the name for the skipped file list (e.g., skipped_files.txt): ").strip()
    
    try:
        max_workers = int(input("Enter the number of parallel workers (e.g., 4 or 8): ").strip())
    except ValueError:
        print("Invalid input for max_workers. Defaulting to 4.")
        max_workers = 4

    run_mafft_on_folder_parallel(
        input_folder=input_folder,
        output_folder=output_folder,
        log_file=log_file,
        skipped_file=skipped_file,
        max_workers=max_workers
    )