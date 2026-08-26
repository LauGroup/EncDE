# License: https://github.com/LauGroup/EncDE/blob/main/LICENSE
# Author: Rezwan Siddiquee; rezwan.siddiquee@sydney.edu.au
# School of Chemistry, The University of Sydney


import pandas as pd
import subprocess
import os # Retained for os.cpu_count() if it were used, or other non-path operations.
from pathlib import Path # Added for path management
import glob # Can be removed if not used elsewhere, current script replaces its usage
from openpyxl import load_workbook
from openpyxl.styles import Font

# --- User Defined Parameters ---
INPUT_DIRECTORY = Path('Demo/')  # Directory containing the Excel files
OUTPUT_DIRECTORY = Path('Demo/') # Folder to save the processed Excel files
CLUSTALO_EXECUTABLE = Path('Scripts/clustalo windows/clustalo.exe') # Path to ClustalO executable

# --- Helper Functions ---

def run_clustalo(sequences, output_fasta_path: Path, iterations=2, full=False, threads=16, force=True):
    # Write sequences to a temporary file for ClustalO input
    temp_input_fasta_path = Path('temp_sequences.fasta') # Temporary file in CWD
    with open(temp_input_fasta_path, 'w') as f:
        for i, seq in enumerate(sequences):
            if seq:  # Ensure the sequence is not empty
                f.write(f">Seq_{i+1}\n{seq}\n")

    print(f"Sequences written to {temp_input_fasta_path}")

    # Construct the Clustal Omega command
    # Ensure paths are strings for subprocess command
    clustalo_cmd_parts = [
        str(CLUSTALO_EXECUTABLE),
        "-i", str(temp_input_fasta_path),
        "-o", str(output_fasta_path),
        "--outfmt=fa",
        f"--iterations={iterations}",
        f"--threads={threads}"
    ]
    if full:
        clustalo_cmd_parts.append("--full")
    if force:
        clustalo_cmd_parts.append("--force")

    clustalo_cmd_str = " ".join(clustalo_cmd_parts) # More robust way to build command string if needed, or pass list directly
    print(f"Running Clustal Omega with command: {clustalo_cmd_str}")

    # Run ClustalO to align the sequences
    # Using shell=False and passing command as a list is generally safer
    result = subprocess.run(clustalo_cmd_parts, capture_output=True, text=True, check=False)

    # Print out the stdout and stderr to help debug
    print("ClustalO stdout:", result.stdout)
    print("ClustalO stderr:", result.stderr)

    if result.returncode != 0:
        # Clean up input fasta before raising error
        if temp_input_fasta_path.exists():
            temp_input_fasta_path.unlink()
        raise Exception(f"ClustalO command failed with exit status {result.returncode}\nstderr: {result.stderr}\nstdout: {result.stdout}")

    # Read the aligned sequences from the output file
    aligned_sequences = []
    current_sequence = ""
    if output_fasta_path.exists():
        with open(output_fasta_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    if current_sequence:
                        aligned_sequences.append(current_sequence)
                    current_sequence = ""
                else:
                    current_sequence += line.strip()
            if current_sequence: # Append the last sequence
                aligned_sequences.append(current_sequence)
    else:
        print(f"Error: ClustalO output file {output_fasta_path} not found.")
        # Clean up input fasta before raising error
        if temp_input_fasta_path.exists():
            temp_input_fasta_path.unlink()
        raise FileNotFoundError(f"ClustalO output file {output_fasta_path} not found after execution.")


    # Clean up temporary files
    if temp_input_fasta_path.exists():
        temp_input_fasta_path.unlink()
        print(f"Removed temporary input file: {temp_input_fasta_path}")
    if output_fasta_path.exists(): # This output_fasta_path is also temporary in this context
        output_fasta_path.unlink()
        print(f"Removed temporary aligned fasta file: {output_fasta_path}")

    return aligned_sequences


def apply_courier_font_to_columns(sheet, columns):
    courier_font = Font(name='Courier New')
    for col_letter in columns: # Changed 'col' to 'col_letter' for clarity if they are 'B', 'C'
        for cell in sheet[col_letter]:
            if cell.value is not None: # Ensure cell has a value before trying to set font
                cell.font = courier_font

def process_excel_file(input_file_path: Path, output_dir_path: Path, iterations=2, full=False, threads=16, force=True):
    print(f"\n--- Starting processing for Excel file: {input_file_path.name} ---")
    output_excel_path = output_dir_path / input_file_path.name
    # Temporary path for ClustalO's aligned fasta output for this specific Excel file
    temp_aligned_fasta_for_excel = output_dir_path / f"{input_file_path.stem}_aligned_temp.fasta"


    try:
        df = pd.read_excel(input_file_path, sheet_name='Sheet 1')

        # Assuming sequences are in the second column (index 1)
        # Clean up sequences: Remove stop codons (*) and trim whitespace
        if len(df.columns) > 1:
            df['Sequence_Cleaned'] = df.iloc[:, 1].astype(str).str.replace('*', '', regex=False).str.strip()
            sequences = df['Sequence_Cleaned'].tolist()
        else:
            print(f"Warning: File '{input_file_path.name}' has less than 2 columns. Cannot extract sequences. Skipping.")
            return

        # Filter out empty sequences before alignment
        non_empty_sequences_for_alignment = [seq for seq in sequences if seq]
        # Keep track of original indices of non-empty sequences if needed, or map results back carefully.
        # For simplicity, we assume ClustalO returns them in order.

        if not non_empty_sequences_for_alignment:
            print(f"No valid sequences found in '{input_file_path.name}' after cleaning. Skipping alignment.")
            # Save the original file to output if desired, or just skip
            # For now, we skip creating an output file if no sequences.
            return

        print(f"Found {len(non_empty_sequences_for_alignment)} non-empty sequences to align from {input_file_path.name}.")
        aligned_sequences = run_clustalo(non_empty_sequences_for_alignment, temp_aligned_fasta_for_excel, iterations, full, threads, force)

        if 'Aligned_Sequence' in df.columns:
            df.drop(columns=['Aligned_Sequence'], inplace=True)

        # Map aligned sequences back. This assumes ClustalO returns aligned sequences
        # in the same order as the non-empty input sequences.
        # Create a new series for aligned sequences, aligning with original DataFrame rows.
        aligned_series = pd.Series([None] * len(df), index=df.index, dtype=object)
        current_aligned_idx = 0
        for i, original_seq_cleaned in enumerate(df['Sequence_Cleaned']):
            if original_seq_cleaned and current_aligned_idx < len(aligned_sequences):
                aligned_series.iloc[i] = aligned_sequences[current_aligned_idx]
                current_aligned_idx += 1

        if len(aligned_sequences) != len(non_empty_sequences_for_alignment):
             print(f"Warning for {input_file_path.name}: Number of aligned sequences ({len(aligned_sequences)}) "
                   f"does not match the number of non-empty input sequences ({len(non_empty_sequences_for_alignment)}). "
                   f"Resulting alignment column might be incomplete or misaligned.")
        
        # Insert the aligned sequences. Assuming original sequences were in col B (index 1),
        # cleaned were temporary, new aligned go to col C (index 2).
        df.insert(2, 'Aligned_Sequence', aligned_series)
        
        # Drop the temporary cleaned sequence column if you added one and don't want it in output
        if 'Sequence_Cleaned' in df.columns:
            df.drop(columns=['Sequence_Cleaned'], inplace=True)


        print(f"Saving processed data to: {output_excel_path}")
        df.to_excel(output_excel_path, sheet_name='Sheet 1', index=False, engine='openpyxl')

        print(f"Applying font changes to: {output_excel_path}")
        wb = load_workbook(output_excel_path)
        sheet = wb['Sheet 1']

        # Determine column letters. If original sequences are in col B (2nd col),
        # and aligned are inserted at index 2 (3rd col), these are B and C.
        # Column B (index 1 in 0-indexed df, openpyxl is 1-indexed)
        # Column C (index 2 in 0-indexed df, openpyxl is 1-indexed)
        # This depends on your actual column structure. The original script hardcoded 'B', 'C'.
        # Let's assume original sequences are in column 2 (letter 'B') and aligned in column 3 (letter 'C')
        columns_to_format = []
        if len(df.columns) > 1: # Original sequence column
            columns_to_format.append(chr(ord('A') + 1)) # Second column (e.g., 'B')
        if 'Aligned_Sequence' in df.columns: # Aligned sequence column
             # Get the column index of 'Aligned_Sequence'
            aligned_col_idx_0based = df.columns.get_loc('Aligned_Sequence')
            columns_to_format.append(chr(ord('A') + aligned_col_idx_0based))

        if columns_to_format:
            apply_courier_font_to_columns(sheet, columns_to_format)
        else:
            print(f"No columns identified for font formatting in {output_excel_path.name}")


        wb.save(output_excel_path)
        print(f"Successfully processed and saved: {output_excel_path.name}")

    except Exception as e:
        print(f"--- ERROR processing file {input_file_path.name}: {e} ---")
        import traceback
        traceback.print_exc()


def process_all_excel_files_in_folder(input_dir: Path, output_dir: Path, iterations=2, full=False, threads=16, force=True):
    print(f"Starting Excel processing...")
    print(f"Input Directory: {input_dir.resolve()}")
    print(f"Output Directory: {output_dir.resolve()}")

    if not input_dir.is_dir():
        print(f"Error: Input directory '{input_dir.resolve()}' does not exist.")
        return

    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        print(f"Ensured output directory exists: {output_dir.resolve()}")
    except OSError as e:
        print(f"Error: Could not create output directory '{output_dir.resolve()}': {e}")
        return

    excel_files = list(input_dir.glob('*.xlsx'))

    if not excel_files:
        print(f"No '.xlsx' files found in '{input_dir.resolve()}'.")
        return

    print(f"Found {len(excel_files)} Excel file(s) to process.")
    processed_count = 0
    error_count = 0

    for file_path_obj in excel_files:
        try:
            process_excel_file(file_path_obj, output_dir, iterations, full, threads, force)
            processed_count += 1
        except Exception as e_main: # Catch errors from process_excel_file if re-raised
            print(f"--- CRITICAL FAILURE processing file '{file_path_obj.name}'. Error: {e_main} ---")
            error_count +=1

    print("\n--- Finished processing all files ---")
    print(f"Successfully processed files: {processed_count}")
    if error_count > 0:
        print(f"Failed to process files: {error_count}")


# --- Main Execution ---
if __name__ == "__main__":
    # Example usage:
    # Parameters for ClustalO can be adjusted here if not passed directly or changed from defaults
    clustalo_iterations = 2
    clustalo_full_distance_matrix = True # Corresponds to --full
    clustalo_threads = 8 # Example: use 8 threads
    clustalo_force_overwrite = True # Corresponds to --force

    # Check if ClustalO executable exists
    if not CLUSTALO_EXECUTABLE.is_file():
        print(f"Error: ClustalO executable not found at {CLUSTALO_EXECUTABLE}")
        print("Please check the CLUSTALO_EXECUTABLE path in the script.")
        exit(1)
    else:
        print(f"Using ClustalO executable: {CLUSTALO_EXECUTABLE}")

    process_all_excel_files_in_folder(
        INPUT_DIRECTORY,
        OUTPUT_DIRECTORY,
        iterations=clustalo_iterations,
        full=clustalo_full_distance_matrix,
        threads=clustalo_threads,
        force=clustalo_force_overwrite
    )