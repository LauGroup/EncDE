# License: https://github.com/LauGroup/EncDE/blob/main/LICENSE
# Author: Rezwan Siddiquee; rezwan.siddiquee@sydney.edu.au
# School of Chemistry, The University of Sydney

import pandas as pd
import glob
import os
import numpy as np
from pathlib import Path

# --- User Defined Parameters ---
INPUT_DIRECTORY = 'Demo/Fastq' # Directory containing the Excel files ('.' for current dir)
OUTPUT_DIRECTORY = 'Demo/' # Folder to save the generated plots

# Define the exact start DNA sequence to search for
START_DNA_SEQUENCE = 'TCCCCGGAGCTGTAC'
# Define the exact end DNA sequence to search for
END_DNA_SEQUENCE = 'ACTGCTGGGGTTTTTCAG'
# Define the expected Wild Type protein sequence
WT_PROTEIN_SEQ = 'SPELYSLLHRVHKDTNVLEIEHVRELITAGVFQ' 
# Define whether to filter out sequences longer than WT_PROTEIN_SEQ
FILTER_LONGER_THAN_WT = True  # Set to True to filter, False to include. On by default.
# Define whether to filter out sequences containing stop codons
FILTER_STOP_CODONS = False     # Set to True to filter, False to include. Off by default.

# --- NEW ---
# Define the exact contaminating peptide motif to filter out
CONTAMINATING_MOTIF = 'GSGGSG'
# Define whether to filter out sequences containing the contaminating motif
FILTER_CONTAMINATING_MOTIF = False # Set to True to filter out, False to include.
# --- End of User Defined Parameters ---

#-------------------------------------------------------------------------------
# Translation Dictionary
translation = {
    "TTT": "F", "TCT": "S", "TAT": "Y", "TGT": "C",
    "TTC": "F", "TCC": "S", "TAC": "Y", "TGC": "C",
    "TTA": "L", "TCA": "S", "TAA": "*", "TGA": "*", # Stop codons
    "TTG": "L", "TCG": "S", "TAG": "*", "TGG": "W", # Stop codons

    "CTT": "L", "CCT": "P", "CAT": "H", "CGT": "R",
    "CTC": "L", "CCC": "P", "CAC": "H", "CGC": "R",
    "CTA": "L", "CCA": "P", "CAA": "Q", "CGA": "R",
    "CTG": "L", "CCG": "P", "CAG": "Q", "CGG": "R",

    "ATT": "I", "ACT": "T", "AAT": "N", "AGT": "S",
    "ATC": "I", "ACC": "T", "AAC": "N", "AGC": "S",
    "ATA": "I", "ACA": "T", "AAA": "K", "AGA": "R",
    "ATG": "M", "ACG": "T", "AAG": "K", "AGG": "R",

    "GTT": "V", "GCT": "A", "GAT": "D", "GGT": "G",
    "GTC": "V", "GCC": "A", "GAC": "D", "GGC": "G",
    "GTA": "V", "GCA": "A", "GAA": "E", "GGA": "G",
    "GTG": "V", "GCG": "A", "GAG": "E", "GGG": "G"
}

def translate_dna(dna_sequence):
    """Translates a DNA sequence into a peptide sequence."""
    peptide = ''
    for i in range(0, len(dna_sequence) - (len(dna_sequence) % 3), 3):
        codon = dna_sequence[i:i+3]
        # Handle codons with 'N' - translate as 'X'
        if 'N' in codon:
             peptide += 'X'
        else:
            # Append amino acid or stop codon ('*')
            peptide += translation.get(codon, '?') # Use '?' for unknown codons 
    return peptide

#-------------------------------------------------------------------------------
# Read the data file (.fastq extension) and return a list of peptide sequences and their frequencies
def SortedPeptideSequencesList(fastqFileLocation, start_dna_seq, end_dna_seq, wt_protein_seq):
    PeptideSequences = {}
    
    try:
        with open(fastqFileLocation, 'r') as RawDataFile:
            Lines = RawDataFile.readlines()
    except FileNotFoundError:
        print(f"Error: File not found - {fastqFileLocation}")
        return []
    except Exception as e:
        print(f"Error reading file {fastqFileLocation}: {e}")
        return []

    processed_sequences = 0
    found_peptides = 0
    skipped_stop_codon = 0
    skipped_longer_than_wt = 0
    skipped_contaminating_motif = 0 # --- NEW ---

    len_wt_protein_seq = len(wt_protein_seq) 

    for i in range(1, len(Lines), 4):
        dna_line = Lines[i].strip()
        processed_sequences += 1
        
        start_index = dna_line.find(start_dna_seq)
        end_index = dna_line.find(end_dna_seq)

        if start_index != -1 and end_index != -1 and start_index < end_index:
            extracted_dna = dna_line[start_index : end_index + len(end_dna_seq)]
            peptide = translate_dna(extracted_dna)
            
            # --- Conditional stop codon check ---
            if FILTER_STOP_CODONS and '*' in peptide:
                skipped_stop_codon += 1
                continue 

            # --- Conditional length check ---
            if FILTER_LONGER_THAN_WT and len(peptide) > len_wt_protein_seq:
                skipped_longer_than_wt += 1
                continue

            # --- NEW: Conditional contaminating motif check ---
            if FILTER_CONTAMINATING_MOTIF and CONTAMINATING_MOTIF in peptide:
                skipped_contaminating_motif += 1
                continue
            # --- END NEW ---

            if peptide: 
                found_peptides += 1
                if peptide not in PeptideSequences:
                    PeptideSequences[peptide] = 1
                else:
                    PeptideSequences[peptide] += 1

    print(f"Processed {processed_sequences} sequences from {fastqFileLocation}.")
    print(f"Found {found_peptides} peptides (that passed all active filters) between specified DNA markers.")
    if FILTER_STOP_CODONS:
        print(f"Skipped {skipped_stop_codon} sequences due to internal stop codons (filter was ON).")
    else:
        # If filter is OFF, skipped_stop_codon counter remains 0. 
        # Sequences with stop codons are processed unless filtered by other means.
        print(f"Stop codon filter was OFF; sequences with stop codons were included (unless filtered by other criteria).")

    if FILTER_LONGER_THAN_WT: 
        print(f"Skipped {skipped_longer_than_wt} sequences for being longer than WT protein sequence (length: {len_wt_protein_seq}, filter was ON).")
    else:
        print(f"Filter for sequences longer than WT was OFF; such sequences were included (unless filtered by other criteria).")
        
    # --- NEW ---
    if FILTER_CONTAMINATING_MOTIF:
        print(f"Skipped {skipped_contaminating_motif} sequences containing the motif '{CONTAMINATING_MOTIF}' (filter was ON).")
    else:
        print(f"Filter for contaminating motif ('{CONTAMINATING_MOTIF}') was OFF; such sequences were included (unless filtered by other criteria).")
    # --- END NEW ---


    PeptideSequencesList = []
    wt_count = PeptideSequences.get(wt_protein_seq, 0) 
    PeptideSequencesList.append(('WT', wt_protein_seq, wt_count))

    for key, value in PeptideSequences.items():
        if key != wt_protein_seq:
            PeptideSequencesList.append((None, key, value))

    SortedPeptideSequences = sorted(PeptideSequencesList[1:], key=lambda x: x[2], reverse=True)
    FinalSortedList = [PeptideSequencesList[0]] + SortedPeptideSequences
    
    return FinalSortedList

#-------------------------------------------------------------------------------
# Create a .XLSX file with a single sheet combining all the data
def SelectionRoundSortedSequenceListGenerator(SortedPeptideSequencesFileName, fastqFileLocation, start_dna_seq, end_dna_seq, wt_protein_seq):
    SortedPeptideSequences = SortedPeptideSequencesList(fastqFileLocation, start_dna_seq, end_dna_seq, wt_protein_seq)
    
    # This warning is generally useful regardless of the filter, as WT is handled specially.
    if '*' in wt_protein_seq:
        print(f"Warning: Defined WT sequence '{wt_protein_seq}' contains a stop codon ('*'). Its counting might be affected by the 'FILTER_STOP_CODONS' setting if it were a read sequence.")
        
    if not SortedPeptideSequences or (len(SortedPeptideSequences) == 1 and SortedPeptideSequences[0][2] == 0 and SortedPeptideSequences[0][1] != wt_protein_seq) : 
         is_only_wt_zero = (len(SortedPeptideSequences) == 1 and SortedPeptideSequences[0][0] == 'WT' and SortedPeptideSequences[0][2] == 0)
         if not SortedPeptideSequences or not is_only_wt_zero:
             print(f"No peptide sequences (that passed all active filters) found or generated for {fastqFileLocation}. Skipping Excel generation.")
             return

    TotalPeptideNumber = sum(Data[2] for Data in SortedPeptideSequences)

    if TotalPeptideNumber == 0 and (not (len(SortedPeptideSequences) == 1 and SortedPeptideSequences[0][0] == 'WT')):
        print(f"Total peptide count (after filters) is zero for {fastqFileLocation}. Skipping Frequency calculation and Excel generation.")
        if any(item[0] == 'WT' for item in SortedPeptideSequences):
             df = pd.DataFrame([item for item in SortedPeptideSequences if item[0]=='WT'], columns=['Name', 'Sequence', 'Reads'])
             df['%Frequency'] = 0.0
        else: 
             df = pd.DataFrame(columns=['Name', 'Sequence', 'Reads', '%Frequency'])
        output_path = f'{SortedPeptideSequencesFileName}.xlsx'
        try:
             df.to_excel(output_path, sheet_name='Sheet 1', index=False)
             print(f"Excel file with WT sequence (zero count) or empty sheet saved to {output_path}")
        except Exception as e:
             print(f"Error writing minimal Excel file {output_path}: {e}")
        return
    elif TotalPeptideNumber == 0 and len(SortedPeptideSequences) == 1 and SortedPeptideSequences[0][0] == 'WT':
        print(f"Only WT sequence found with zero counts (after filters) for {fastqFileLocation}. Generating minimal Excel.")
        df = pd.DataFrame([SortedPeptideSequences[0]], columns=['Name', 'Sequence', 'Reads'])
        df['%Frequency'] = 0.0
        Reads_header = f"Reads (# of seq from total 0 seq)"
        Frequency_header = f"%Frequency (# of seq/total 0 seq)"
        df.rename(columns={'Reads': Reads_header, '%Frequency': Frequency_header}, inplace=True)
        df[f"Reads (# of seq from total 0 seq) after WT filtered out"] = np.nan
        df[f"%Frequency (# of seq/total 0 seq) after WT filtered out"] = np.nan
        df[f"Reads (vs WT count 0)"] = 0
        df[f"%Frequency (vs WT count 0)"] = 0.0 

        output_path = f'{SortedPeptideSequencesFileName}.xlsx'
        try:
            df.to_excel(output_path, sheet_name='Sheet 1', index=False)
            print(f"Excel file with WT sequence (zero count) saved to {output_path}")
        except Exception as e:
            print(f"Error writing minimal Excel file {output_path}: {e}")
        return

    Reads_header = f"Reads (# of seq from total {TotalPeptideNumber} seq)"
    Frequency_header = f"%Frequency (# of seq/total {TotalPeptideNumber} seq)"

    enriched_data = []
    for seq_data in SortedPeptideSequences:
        Frequency_value = (seq_data[2] / TotalPeptideNumber * 100) if TotalPeptideNumber > 0 else 0
        enriched_data.append((*seq_data, Frequency_value))

    df = pd.DataFrame(enriched_data, columns=['Name', 'Sequence', Reads_header, Frequency_header])

    non_wt_indices = df[df['Name'] != 'WT'].index
    name_map = {index: i + 1 for i, index in enumerate(non_wt_indices)}
    df['Name'] = df.apply(lambda row: 'WT' if row['Name'] == 'WT' else name_map.get(row.name), axis=1)

    Reads_header_noWT = "N/A"
    Frequency_header_noWT = "N/A"
    Reads_header_overWT = "N/A"
    Frequency_header_overWT = "N/A"

    wt_row = df[df['Name'] == 'WT']
    wt_Reads = wt_row[Reads_header].iloc[0] if not wt_row.empty else 0

    if wt_Reads > 0:
        TotalPeptideNumber_NoWT = df[df['Name'] != 'WT'][Reads_header].sum()
        
        if TotalPeptideNumber_NoWT > 0:
             Reads_header_noWT = f"Reads (# of seq from total {TotalPeptideNumber_NoWT} seq) after WT filtered out"
             Frequency_header_noWT = f"%Frequency (# of seq/total {TotalPeptideNumber_NoWT} seq) after WT filtered out"
             df[Reads_header_noWT] = df.apply(lambda row: np.nan if row['Name'] == 'WT' else row[Reads_header], axis=1)
             df[Frequency_header_noWT] = df.apply(lambda row: np.nan if row['Name'] == 'WT' else (row[Reads_header] / TotalPeptideNumber_NoWT * 100), axis=1)
        else: 
             df[Reads_header_noWT] = np.nan
             df[Frequency_header_noWT] = np.nan
             Reads_header_noWT = f"Reads (# of seq from total 0 seq) after WT filtered out"
             Frequency_header_noWT = f"%Frequency (# of seq/total 0 seq) after WT filtered out"

        Reads_header_overWT = f"Reads (vs WT count {wt_Reads})" 
        Frequency_header_overWT = f"%Frequency (vs WT count {wt_Reads})" 
        df[Reads_header_overWT] = df[Reads_header] 
        df[Frequency_header_overWT] = df.apply(lambda row: 100.0 if row['Name'] == 'WT' else (row[Reads_header] / wt_Reads * 100), axis=1)
        
    else: 
        TotalPeptideNumber_NoWT = df[df['Name'] != 'WT'][Reads_header].sum()
        
        if TotalPeptideNumber_NoWT > 0 :
             Reads_header_noWT = f"Reads (# of seq from total {TotalPeptideNumber_NoWT} seq) after WT filtered out"
             Frequency_header_noWT = f"%Frequency (# of seq/total {TotalPeptideNumber_NoWT} seq) after WT filtered out"
             df[Reads_header_noWT] = df.apply(lambda row: np.nan if row['Name'] == 'WT' else row[Reads_header], axis=1)
             df[Frequency_header_noWT] = df.apply(lambda row: np.nan if row['Name'] == 'WT' else (row[Reads_header] / TotalPeptideNumber_NoWT * 100), axis=1)
        else: 
             df[Reads_header_noWT] = np.nan
             df[Frequency_header_noWT] = np.nan
             Reads_header_noWT = f"Reads (# of seq from total 0 seq) after WT filtered out"
             Frequency_header_noWT = f"%Frequency (# of seq/total 0 seq) after WT filtered out"

        df[Reads_header_overWT] = df[Reads_header] 
        Frequency_header_overWT = f"%Frequency (vs WT count 0)" 
        df[Frequency_header_overWT] = np.nan 
        df.loc[df['Name'] == 'WT', Frequency_header_overWT] = 0.0

    final_columns = [
        'Name',
        'Sequence',
        Reads_header,
        Frequency_header,
        Reads_header_noWT,
        Frequency_header_noWT,
        Reads_header_overWT,
        Frequency_header_overWT
    ]
    
    for col in final_columns:
        if col not in df.columns and col != "N/A":
            df[col] = np.nan
            
    final_columns_filtered = [col for col in final_columns if col != "N/A"]
    df = df[final_columns_filtered]

    output_path = f'{SortedPeptideSequencesFileName}.xlsx'
    try:
        df.to_excel(output_path, sheet_name='Sheet 1', index=False)
        print(f"Excel file with combined sheet saved to {output_path}")
    except Exception as e:
        print(f"Error writing Excel file {output_path}: {e}")

#-------------------------------------------------------------------------------
# Process all fastq files in the defined directory
if __name__ == "__main__":
    print(f"Starting analysis with:")
    print(f"  Start DNA Sequence: {START_DNA_SEQUENCE}")
    print(f"  End DNA Sequence:   {END_DNA_SEQUENCE}")
    print(f"  WT Protein Sequence: {WT_PROTEIN_SEQ}")
    print(f"  Filter sequences longer than WT: {'On' if FILTER_LONGER_THAN_WT else 'Off'}")
    print(f"  Filter sequences with stop codons: {'On' if FILTER_STOP_CODONS else 'Off (stop codons included)'}")
    # --- NEW ---
    print(f"  Filter sequences with motif '{CONTAMINATING_MOTIF}': {'On' if FILTER_CONTAMINATING_MOTIF else 'Off (motif-containing sequences included)'}")
    # --- END NEW ---
    
    script_dir = Path(INPUT_DIRECTORY).resolve()
    print(f"Searching for .fastq files in: {script_dir}")
    
    output_dir = Path(OUTPUT_DIRECTORY)
    output_dir.mkdir(parents=True, exist_ok=True)

    fastq_files = glob.glob(os.path.join(script_dir, "*.fastq"))

    if not fastq_files:
        print("No .fastq files found in the script directory.")
    else:
        print(f"Found {len(fastq_files)} .fastq files: {fastq_files}")
        for fastq_file_path_str in fastq_files:
            fastq_file_path = Path(fastq_file_path_str)
            print(f"\nProcessing file: {fastq_file_path.name}")
            output_file_name = fastq_file_path.stem
            excel_file_path_prefix = output_dir/output_file_name
            

            SelectionRoundSortedSequenceListGenerator(str(excel_file_path_prefix), str(fastq_file_path), START_DNA_SEQUENCE, END_DNA_SEQUENCE, WT_PROTEIN_SEQ)

    print("\nAnalysis complete.")