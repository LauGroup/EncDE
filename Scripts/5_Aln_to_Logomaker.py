# License: https://github.com/LauGroup/EncDE/blob/main/LICENSE
# Author: Rezwan Siddiquee; rezwan.siddiquee@sydney.edu.au
# School of Chemistry, The University of Sydney


import os
import pandas as pd
import logomaker
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np # For information content calculation if needed

# --- User Defined Parameters ---

# --- Directory Settings ---
INPUT_DIRECTORY = Path('Logos/')  # Directory containing the Excel files
OUTPUT_DIRECTORY = Path('Logos/') # Folder to save the generated logo images

# --- Logo Data & Weighting ---
SEQUENCE_COLUMN_INDEX = 2
WT_SEQUENCE_COLUMN_INDEX = 2
FREQUENCY_COLUMN_INDEX = 5
USE_FREQUENCY_WEIGHTING = False

# --- Logo Sequence Slicing (Multiple Ranges) ---
LOGO_SEQUENCE_RANGES = [(1, 33)] # Example: [(1, 5), (7, 10), (12, 15)] (1-based, inclusive) (1,33 is full sequence)
                          # If empty or None, the full sequence will be used.

# --- Y-axis Configuration ---
Y_AXIS_TYPE = 'bits'
SHOW_Y_AXIS = True
Y_AXIS_LABEL_BITS = 'Frequency'
Y_AXIS_LABEL_PROBABILITY = 'Probability'
Y_AXIS_FONT_NAME = None
Y_AXIS_FONT_SIZE = 5

# --- X-axis Configuration ---
# Bottom X-axis (WT Sequence)
SHOW_BOTTOM_X_AXIS = True
BOTTOM_X_AXIS_FONT_NAME = None
BOTTOM_X_AXIS_FONT_SIZE = 10

# Top X-axis (Continuous Numbering)
SHOW_TOP_X_AXIS = True
# Starting number for the first position of the displayed logo on the top X-axis.
TOP_X_AXIS_NUMBERING_START = 185 # Default: 1
# Interval for top x-axis ticks (e.g., every Nth position of the *final concatenated logo*)
TOP_X_AXIS_TICK_INTERVAL = 5
TOP_X_AXIS_FONT_NAME = None
TOP_X_AXIS_FONT_SIZE = 10

# --- Plot Aesthetics ---
PLOT_WIDTH = 5
PLOT_HEIGHT = 1
PLOT_TITLE_PREFIX = "Sequence Logo for "
TITLE_FONT_SIZE = 12
PLOT_TOP_MARGIN = 0.85
LOGO_FONT_NAME = 'Arial'
LOGO_COLOR_SCHEME = 'weblogo_protein'
# --- End User Defined Parameters ---

def get_concatenated_sliced_data(full_sequence, ranges_1_based_inclusive, original_length_for_validation):
    """
    Slices a sequence based on multiple 1-based inclusive ranges and concatenates them.
    Also returns a map of new indices to original 1-based indices (though this map is not used for axis numbering in this version).
    """
    if not ranges_1_based_inclusive: # If no ranges, use full sequence
        sliced_parts = [full_sequence]
        # original_indices_map for full sequence (1 to length)
        original_indices_map = list(range(1, len(full_sequence) + 1)) if full_sequence else []
        return "".join(sliced_parts), original_indices_map

    sliced_parts = []
    original_indices_map = []

    for r_start, r_end in ranges_1_based_inclusive:
        if not (isinstance(r_start, int) and isinstance(r_end, int) and r_start >= 1 and r_end >= r_start):
            print(f"Warning: Invalid range ({r_start}, {r_end}). Skipping this range.")
            continue
        
        py_start = r_start - 1
        py_end = r_end

        if py_start >= original_length_for_validation:
            continue

        actual_py_end = min(py_end, original_length_for_validation)
        if py_start >= actual_py_end:
            continue

        segment = full_sequence[py_start:actual_py_end]
        sliced_parts.append(segment)

        for i in range(r_start, r_start + len(segment)):
            original_indices_map.append(i)

    return "".join(sliced_parts), original_indices_map


def generate_logo_from_excel(file_path_obj, config):
    try:
        xls = pd.ExcelFile(file_path_obj)
        df = xls.parse(xls.sheet_names[0])
    except Exception as e:
        print(f"Error reading Excel file {file_path_obj.name}: {e}")
        return

    try:
        wt_sequence_full = str(df.iloc[0, config.WT_SEQUENCE_COLUMN_INDEX]) if config.WT_SEQUENCE_COLUMN_INDEX < len(df.columns) and not df.empty else ""
        sequences_full_list = df.iloc[:, config.SEQUENCE_COLUMN_INDEX].astype(str).tolist() if config.SEQUENCE_COLUMN_INDEX < len(df.columns) else []
        
        if not sequences_full_list:
             print(f"Error: Sequence column {config.SEQUENCE_COLUMN_INDEX} not found or empty in {file_path_obj.name}. Skipping.")
             return

        counts = None
        current_use_weighting = config.USE_FREQUENCY_WEIGHTING
        if current_use_weighting:
            if config.FREQUENCY_COLUMN_INDEX < len(df.columns):
                counts = df.iloc[:, config.FREQUENCY_COLUMN_INDEX].fillna(0).astype(int).tolist()
            else:
                print(f"Warning: Freq column {config.FREQUENCY_COLUMN_INDEX} not found in {file_path_obj.name}. Weighting disabled.")
                current_use_weighting = False
    except Exception as e:
        print(f"Error extracting data columns from {file_path_obj.name}: {e}")
        return

    # --- Process Slicing based on LOGO_SEQUENCE_RANGES ---
    processed_sequences = []
    # The original_indices_map is calculated but not used for top X-axis numbering in this version
    wt_sequence_sliced, _ = get_concatenated_sliced_data(wt_sequence_full, config.LOGO_SEQUENCE_RANGES, len(wt_sequence_full))

    for seq_full in sequences_full_list:
        sliced_seq, _ = get_concatenated_sliced_data(seq_full, config.LOGO_SEQUENCE_RANGES, len(seq_full))
        processed_sequences.append(sliced_seq)

    if not any(processed_sequences):
        print(f"Warning: All sequences are empty after slicing for {file_path_obj.name}. Skipping logo.")
        return
    
    sequence_list_for_matrix = []
    if current_use_weighting and counts and len(counts) == len(processed_sequences):
        for seq, count_val in zip(processed_sequences, counts):
            if seq: sequence_list_for_matrix.extend([seq] * int(count_val))
    else:
        for seq in processed_sequences:
            if seq: sequence_list_for_matrix.append(seq)

    if not sequence_list_for_matrix:
        print(f"No valid sequences to process for {file_path_obj.name} after weighting/slicing. Skipping.")
        return

    matrix_len = 0
    if sequence_list_for_matrix:
         matrix_len = len(sequence_list_for_matrix[0])
         if matrix_len == 0:
             print(f"Sequences for matrix in {file_path_obj.name} have zero length. Skipping.")
             return
         if not all(len(s) == matrix_len for s in sequence_list_for_matrix):
             print(f"Error: Not all sequences for matrix have same length ({matrix_len}) in {file_path_obj.name}. Skipping.")
             return
    else:
        print(f"No sequences for matrix in {file_path_obj.name}. Skipping.")
        return

    counts_matrix = logomaker.alignment_to_matrix(sequence_list_for_matrix, to_type='counts')
    y_label = config.Y_AXIS_LABEL_BITS
    if config.Y_AXIS_TYPE == 'probability':
        matrix_to_plot = logomaker.transform_matrix(counts_matrix, from_type='counts', to_type='probability')
        y_label = config.Y_AXIS_LABEL_PROBABILITY
    else: # Default or 'bits'
        matrix_to_plot = counts_matrix
        if config.Y_AXIS_TYPE != 'bits':
             print(f"Warning: Unknown Y_AXIS_TYPE '{config.Y_AXIS_TYPE}'. Defaulting to 'bits'.")

    fig, ax = plt.subplots(figsize=(config.PLOT_WIDTH, config.PLOT_HEIGHT))
    fig.subplots_adjust(top=config.PLOT_TOP_MARGIN)

    try:
        logo = logomaker.Logo(matrix_to_plot,
                              ax=ax,
                              color_scheme=config.LOGO_COLOR_SCHEME,
                              font_name=config.LOGO_FONT_NAME)
    except Exception as e:
        print(f"Error creating Logomaker logo for {file_path_obj.name}: {e}")
        plt.close(fig)
        return

    if config.SHOW_Y_AXIS:
        ax.set_ylabel(y_label, fontname=config.Y_AXIS_FONT_NAME, fontsize=config.Y_AXIS_FONT_SIZE)
        ax.tick_params(axis='y', labelsize=config.Y_AXIS_FONT_SIZE, labelrotation=0)
        for ticklabel in ax.get_yticklabels():
            ticklabel.set_fontname(config.Y_AXIS_FONT_NAME)
    else:
        ax.set_ylabel("")
        ax.set_yticks([])

    ax.set_xlabel("")
    if config.SHOW_BOTTOM_X_AXIS and wt_sequence_sliced:
        ax.set_xticks(range(len(wt_sequence_sliced)))
        ax.set_xticklabels(list(wt_sequence_sliced),
                           fontname=config.BOTTOM_X_AXIS_FONT_NAME,
                           fontsize=config.BOTTOM_X_AXIS_FONT_SIZE)
    else:
        ax.set_xticks([])
        ax.set_xticklabels([])

    ax2 = None
    if config.SHOW_TOP_X_AXIS and matrix_len > 0:
        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())

        # Continuous numbering for top X-axis
        tick_positions_top = range(0, matrix_len, config.TOP_X_AXIS_TICK_INTERVAL)
        # Labels start from TOP_X_AXIS_NUMBERING_START and increment by 1 for each position in the logo
        number_labels_top = [str(config.TOP_X_AXIS_NUMBERING_START + tick_pos) for tick_pos in tick_positions_top]
        
        ax2.set_xticks(tick_positions_top)
        ax2.set_xticklabels(number_labels_top,
                            fontname=config.TOP_X_AXIS_FONT_NAME,
                            fontsize=config.TOP_X_AXIS_FONT_SIZE)
        ax2.set_xlabel("")
    elif ax2: # Should not be created if not shown, but as a safeguard
        ax2.set_visible(False)

    plot_title_text = f"{config.PLOT_TITLE_PREFIX}{file_path_obj.stem}"
    ax.set_title(plot_title_text, fontsize=config.TITLE_FONT_SIZE)

    if not config.OUTPUT_DIRECTORY.exists():
        config.OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    base_output_path = config.OUTPUT_DIRECTORY / file_path_obj.stem

    try:
        plt.savefig(f"{base_output_path}.pdf", format='pdf', bbox_inches='tight')
        print(f"Saved {base_output_path}.pdf")
    except Exception as e:
        print(f"Error saving logo for {file_path_obj.name}: {e}")
    finally:
        plt.close(fig)

# --- Main Script Execution ---
if __name__ == "__main__":
    config_params = {
        "INPUT_DIRECTORY": INPUT_DIRECTORY, "OUTPUT_DIRECTORY": OUTPUT_DIRECTORY,
        "SEQUENCE_COLUMN_INDEX": SEQUENCE_COLUMN_INDEX, "WT_SEQUENCE_COLUMN_INDEX": WT_SEQUENCE_COLUMN_INDEX,
        "FREQUENCY_COLUMN_INDEX": FREQUENCY_COLUMN_INDEX, "USE_FREQUENCY_WEIGHTING": USE_FREQUENCY_WEIGHTING,
        "LOGO_SEQUENCE_RANGES": LOGO_SEQUENCE_RANGES,
        "Y_AXIS_TYPE": Y_AXIS_TYPE, "SHOW_Y_AXIS": SHOW_Y_AXIS,
        "Y_AXIS_LABEL_BITS": Y_AXIS_LABEL_BITS, "Y_AXIS_LABEL_PROBABILITY": Y_AXIS_LABEL_PROBABILITY,
        "Y_AXIS_FONT_NAME": Y_AXIS_FONT_NAME, "Y_AXIS_FONT_SIZE": Y_AXIS_FONT_SIZE,
        "SHOW_BOTTOM_X_AXIS": SHOW_BOTTOM_X_AXIS,
        "BOTTOM_X_AXIS_FONT_NAME": BOTTOM_X_AXIS_FONT_NAME, "BOTTOM_X_AXIS_FONT_SIZE": BOTTOM_X_AXIS_FONT_SIZE,
        "SHOW_TOP_X_AXIS": SHOW_TOP_X_AXIS,
        "TOP_X_AXIS_NUMBERING_START": TOP_X_AXIS_NUMBERING_START, # Added back
        "TOP_X_AXIS_TICK_INTERVAL": TOP_X_AXIS_TICK_INTERVAL,
        "TOP_X_AXIS_FONT_NAME": TOP_X_AXIS_FONT_NAME, "TOP_X_AXIS_FONT_SIZE": TOP_X_AXIS_FONT_SIZE,
        "PLOT_WIDTH": PLOT_WIDTH, "PLOT_HEIGHT": PLOT_HEIGHT,
        "PLOT_TITLE_PREFIX": PLOT_TITLE_PREFIX, "TITLE_FONT_SIZE": TITLE_FONT_SIZE,
        "PLOT_TOP_MARGIN": PLOT_TOP_MARGIN, "LOGO_FONT_NAME": LOGO_FONT_NAME,
        "LOGO_COLOR_SCHEME": LOGO_COLOR_SCHEME,
    }

    class Config:
        def __init__(self, entries): self.__dict__.update(entries)
    current_config = Config(config_params)

    if not current_config.INPUT_DIRECTORY.exists() or not current_config.INPUT_DIRECTORY.is_dir():
        print(f"Error: Input directory '{current_config.INPUT_DIRECTORY}' not found.")
    else:
        print(f"Scanning for Excel files in: {current_config.INPUT_DIRECTORY.resolve()}")
        excel_files_found = list(current_config.INPUT_DIRECTORY.glob("*.xlsx"))
        if not excel_files_found:
            print(f"No Excel files (.xlsx) found in '{current_config.INPUT_DIRECTORY}'.")
        else:
            print(f"Found {len(excel_files_found)} Excel file(s). Processing...")
            for file_path_obj in excel_files_found:
                print(f"\nProcessing file: {file_path_obj.name}")
                generate_logo_from_excel(file_path_obj, current_config)
            print("\nAll files processed.")