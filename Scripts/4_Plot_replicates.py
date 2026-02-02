# License: https://github.com/LauGroup/EncDE/blob/main/LICENSE
# Author: Rezwan Siddiquee; rezwan.siddiquee@sydney.edu.au
# School of Chemistry, The University of Sydney

import pandas as pd
import os
import glob
import time
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.cm as cm
import re
import itertools
import math
import numpy as np
from pathlib import Path
import warnings
from scipy.stats import linregress # Added for R-squared calculation

# ==============================================================================
#                            CONFIGURATION SETTINGS
# ==============================================================================
# --- General Parameters ---
INPUT_DIRECTORY = 'Demo/'  # Directory containing the Excel files ('.' for current dir)
OUTPUT_DIRECTORY = 'Demo/'  # Folder to save the generated plots and summary Excel

# --- Plotting & Analysis Switches ---
PLOT_RAW_Frequency = False       # Generate Raw Frequency plots
PLOT_LOG2_Frequency = True      # Generate Log2 Raw Frequency plots
PLOT_FOLD_CHANGE = False          # Generate Fold Change plots
PLOT_LOG2_FOLD_CHANGE = True     # Generate Log2 Fold Change plots
CREATE_EXCEL_OUTPUT_SUMMARY = False # Create the Excel summary file

# --- Fold Change Parameters ---
# This new logic uses Reads Per Million (RPM). A baseline file (e.g., "1_Input.xlsx") is required.
# FC = RPM_selection / RPM_input.
USE_FILTER_THRESHOLD = True        # If True, applies a different FC calculation for sequences with low reads in the selection sample.
FILTER_READ_THRESHOLD = 70         # Read count threshold for the filter.
                                   # If USE_FILTER_THRESHOLD is True, sequences in the selection with reads < this value will have their
                                   # FC calculated as: FC = RPM_selection / Filter_RPM.
                                   # "Filter_RPM" is the highest RPM found among sequences that have *exactly* FILTER_READ_THRESHOLD reads in the same selection sample.


# --- Data Filtering ---
INCLUDE_WT_IN_OUTPUT = True      # If False, the specified WT sequence will be excluded from all plots and Excel summaries.
                                   # The WT sequence is defined in the STYLING_RULES section below.
FILTER_BY_DELETION_COUNT = True    # NEW: If True, filters out sequences with a deletion count <= the threshold below.
DELETION_COUNT_THRESHOLD = 0       # NEW: Threshold for deletion filtering. WT is always kept, regardless of this filter.

# --- File Name Parsing ---
# Regex to extract info from filename (e.g., "2A_Selection1_R1.xlsx")
# Group 1: Condition Number (e.g., "2")
# Group 2: Replicate ID (e.g., "A")
# Group 3: Condition Label (e.g., "Selection1_R1")
FILENAME_REGEX = r'^(\d+)([A-Za-z]?)_?(.*?)(\.xlsx)?$'

# --- Plotting Parameters (individual plot aesthetics) ---
PLOT_APPEARANCE_CONFIG = {
    # --- Data Selection for individual plots ---
    # 'axis_type' will be set dynamically based on the above PLOT_ boolean switches
    'plot_top_n': None,  # Set to an integer to plot only top N sequences, None to plot all.
                         # Ranking based on average value of the chosen axis_type for the plotted pair.

    # --- Figure Aesthetics ---
    'fig_size': (3.4, 1.6),  # IF not showing seq in legend use this.
    #'fig_size': (3, 3),  # If showing seq in legend use this.
    'font_family': 'Arial',
    'font_size': 7,
    'dot_size': 10,

    # --- Axis Limits ---
    #'x_limit': None,  # Set to a tuple (min, max) for manual limits, None for auto.
    #'y_limit': None,  # Set to a tuple (min, max) for manual limits, None for auto.
    'x_limit': (None),  # Set to a tuple (min, max) for manual limits, None for auto.
    'y_limit': (None),  # Set to a tuple (min, max) for manual limits, None for auto.
    # --- Output Formats ---
    'save_formats': ['pdf'],  # File formats to save plots ('pdf', 'svg', 'png')
    'dpi': 300,  # Resolution for rasterized formats (png)

    # --- Top N Sequences in Legend ---
    'SHOW_TOP_SEQUENCES_IN_LEGEND': True, # Toggle to show top N sequences in legend
    'LEGEND_TOP_N_SEQUENCES': 0,          # Number of top sequences to list if above is True
    'TOP_SEQUENCE_MARKER_SIZE_FACTOR': 1, # Factor to increase marker size for top N sequences
    'SHOW_R_SQUARE_IN_LEGEND': True,      # NEW: Toggle to show R-squared value in legend
}

# --- Sequence Styling Rules ---
# Rules are applied in order (index 0 has highest priority).
WT_SEQUENCE = "SPELYSLLHRVHKDTNVLEIEHVRELITAGVFQ" # Define WT sequence here for easy access


YOLO4 = "SPELYSLLHRVHKDYIEHVRELITAGVFQ"

SLUGS8 = "SPELYSLLGSGGSGVRELITAGVFQ"
HOLEBOY13 = "SPELYGSGGSGELITAGVFQ"

PIGS13 = "SPELYPYGSCIELITAGVFQ"

GLASS9 = "SPELYSLLHRGSGGSGLITAGVFQ"
LETTER11 = "SPELYSLLGSGGSGLITAGVFQ"
SLAY13 = "SPELYSLGSGGSGITAGVFQ"
VEGETA10 = "SPELYSLLHRVGSGGSGTAGVFQ"

# BASE_STYLING_RULES define the primary color/style for a dot.
# The highest-priority (first) matching rule will be used.
BASE_STYLING_RULES = [
    {
        'condition': {'has_stop_codon': True},
        'style': {'label': 'Stop Codons', 'color': '#EB1E4E', 'alpha': 1, 'edgecolor': '#000000', 'linewidth': 0, 'zorder': 0}
    },
    {
        'condition': {'is_frameshift': True},
        'style': {'label': 'Frameshifts', 'color': '#EB1E4E', 'alpha': 1, 'edgecolor': '#000000', 'linewidth': 0, 'zorder': 0}
    },
    {
        'condition': {'deletion_range': (0, 4)},
        'style': {'label': '0-4 Deletions', 'color': '#f0f921', 'alpha': 0.7, 'zorder': 1}
    },
    {
        'condition': {'deletion_range': (5, 12)},
        'style': {'label': '5-12 Deletions', 'color': '#ed7953', 'alpha': 0.7, 'zorder': 1}
    },
    {
        'condition': {'deletion_range': (13, 13)},
        'style': {'label': '13 Deletions', 'color': '#9c179e', 'alpha': 0.7, 'zorder': 1}
    },
    {
        'condition': {'deletion_range': (14, 33)},
        'style': {'label': '14+ Deletions', 'color': '#0d0887', 'alpha': 0.7, 'zorder': 1}
    },
    {
        'condition': {}, # Matches all sequences (lowest priority)
        'style': {'label': 'Other', 'color': '#808080', 'alpha': 0.3, 'edgecolor': None, 'linewidth': 0}
    },
]

# OVERLAY_STYLING_RULES are applied *after* the base rule.
# ALL matching rules in this list will be applied, updating the style.
# This allows you to add/override properties like 'edgecolor' or 'zorder'.
OVERLAY_STYLING_RULES = [
    {
        'condition': {'sequence_list': [WT_SEQUENCE]}, # Reference the WT sequence variable
        # This rule *includes* a color, so it will override the base deletion color.
        'style': {'label': 'WT', 'color': '#000000', 'alpha': 1, 'edgecolor': '#000000', 'linewidth': 1, 'zorder': 5}
    },
    {
        'condition': {'sequence_list': []}, # Reference the sequence variable
        # This rule *omits* 'color', so it will keep the base deletion color for the plot
        # but apply the specified alpha, edgecolor, etc.
        'style': {
            'label': 'Selected', 
            'alpha': 1, 
            'edgecolor': '#000000', 
            'linewidth': 1, 
            'zorder': 5,
            'legend_color': "#FFFFFF" # <<< EDIT: Specific color for the legend entry
        }
    },
    
]


# ==============================================================
#                  PREDEFINED SEQUENCE NAMES
# ==============================================================
predefined_names = {
    "SPELYGSGGSGELITAGVFQ": "HoleBoy13",
    "SPELYSLLGSGGSGVRELITAGVFQ": "Slugs8",
    "SPELYSLLHRVHKHTYVLQIEHVRELITAGVFQ": "HYQ",
    "SPELYSLLHRVHKDTNVLSIEHVRELITAGVFQ": "Sieve",
    "SPELYSLLHRVHKHYIQHVRELITAGVFQ": "Hookie4",
    "SPELYSLLHRVHKDYIEHVRELITAGVFQ": "Yolo4",
    "SPELYSLLGSGGSGEHVRELITAGVFQ": "Every6",
    "SPELYGSGGSGVLITAGVFQ": "Volvo13",
    "SPELYSLLHRVHKDTNVLEIEHVRELITAGVFQ": "WT",
    "SPELYSLLHRVHLITAGVFQ": "Holiday13",
    "SPELYTMFFFFELITAGVFQ": "Feli13",
    "SPELYSLLHRVHSVRELITAGVFQ": "Vest9",
    "SPELYRLCSFRELITAGVFQ": "Lyrics13",
    "SPELYGHAYTLELITAGVFQ": "Guy13",
    "SPELYFHMLFPELITAGVFQ": "Fam13",
    "SPELYSLLHRVHKDTNVLSIEHVCMVLLGFF": "Slice2",
    "SPELYHYLTSPELITAGVFQ": "Hotel13",
    "SPELYSLLHRVHWELITAGVFQ": "Weli11",
    "SPELYNIEHVRELITAGVFQ": "Never13",
    "SPELYSLGSGGSGLITAGVFQ": "Legs",
    "SPELYSLLHRVHRITAGVFQ": "Slave",
    "SPELYGSGGSGITAGVFQ": "Gita",
    "SPELYSLLHRVHWDTNVLEIEHVRELITAGVFQ": "Wind",
    "SPELYSLLGSGGSGLITAGVFQ": "Letter11",
    "SPELYSLLGSGGSGITAGVFQ": "Guest",
    "SPELYSLLHRGSGGSGITAGVFQ": "Harvey",
    "SPELYSLLHRVHGTAGVFQ": "Hagrid",
    "SPELYSLGSGGSGRELITAGVFQ": "Rat",
    "SPELYSLLHRVHLTAGVFQ": "Tail",
    "SPELYSLLHRVHRTAGVFQ": "Shell",
    "SPELYGSGGSGLITAGVFQ": "Glide",
    "SPELYSLLGSGGSGRELITAGVFQ": "Greta",
    "SPELYSLIIEHVRELITAGVFQ": "Sly",
    "SPELYGSGGSGRELITAGVFQ": "Grey",
    "SPELYSLLHRVHKGSGGSGLITAGVFQ": "Hike",
    "SPELYSGSGGSGRELITAGVFQ": "Really",
    "SPELYSLLHRVGSGGSGRELITAGVFQ": "Ruby",
    "SPELYSLLHRVHKDNIEHVRELITAGVFQ": "Dan",
    "SPELYSLLHRVHKGSGGSGELITAGVFQ": "Vicky",
    "SPELYSLLGSGGSGEIEHVRELITAGVFQ": "Geisha",
    "SPELYSLLHRVGSGGSGHVRELITAGVFQ": "Heavy",
    "SPELYSLLHRGSGGSGVRELITAGVFQ": "Hugs",
    "SPELYSLLHRVGSGGSGLEIEHVRELITAGVFQ": "Glee",
    "SPELYSGSGGSGNVLEIEHVRELITAGVFQ": "Gwen",
    "SPELYSLLHRVQWELITAGVFQ": "Qweli",
    "SPELYSLLHRVHKDTNVLEIVHVRELITAGVFQ": "Live",
    "SPELYSLLHRVGSGGSGLITAGVFQ": "River",
    "SPELYHYLSSPELITAGVFQ": "Special",
    "SPELYTMFFFFEMITAGVFQ": "Femi",
    "SPELYSLHHPNELITAGVFQ": "Neli",
    "SPELYSLLHRGSGGSGLITAGVFQ": "Glass",
    "SPELYSLLGSGGSGELITAGVFQ": "Glimmer",
    "SPELYSLLHGSGGSGLEIEHVRELITAGVFQ": "Huge",
    "SPELYSLGSGGSGITAGVFQ": "Slay13",
    "SPELYSLLGSGGSGTAGVFQ": "Slim",
    "SPELYPYGSCIELITAGVFQ": "Pigs13",
    "SPELYSLLHRVGSGGSGTAGVFQ": "Vegeta10",
    "SPELYGSGSSGELITAGVFQ": "Ellie",
    "SPELYSLLGCGGSGVRELITAGVFQ": "Give",
    "SPELYSLLHHVHKDTNVLEIEHVRELITAGVFQ": "Hunter",
    "SPELYSLLHHVHKDYIEHVRELITAGVFQ": "Holo",
    "SPELYSLSHRVHKDYIEHVRELITAGVFQ": "Solo",
    "SPELYSLGSGGSGELITAGVFQ": "Elsa",
    "SPELYSLLHGSGGSGITAGVFQ": "Hello",
}
# ==============================================================================
#                        HELPER FUNCTIONS
# ==============================================================================

def parse_filename(filename):
    """Parses filename to extract condition number, replicate ID, and label."""
    match = re.match(FILENAME_REGEX, os.path.basename(filename))
    if match:
        num_str, rep_id, label, _ = match.groups()
        try:
            num = int(num_str)
            rep_id = rep_id if rep_id else None # Keep as None if empty, not ""
            label = label.replace('.xlsx', '').strip('_')
            return num, rep_id, label
        except ValueError:
            return None, None, None
    else:
        print(f"Warning: Could not parse filename: {filename} using regex: {FILENAME_REGEX}")
        return None, None, None

def calculate_log2_value(value):
    """Calculates log2 of a value, handling non-positive values."""
    if value > 0:
        return np.log2(value)
    elif value == 0:
        return -np.inf # Represent log2(0) as negative infinity
    else: # value < 0 or NaN
        return np.nan # Log2 of negative number or NaN is undefined/NaN

# === EDIT START: Updating get_style_for_sequence function ===
def get_style_for_sequence(sequence, x_value, y_value, deletion_count, config):
    """
    Determines the plot style for a sequence based on prioritized rules.
    Returns:
    - final_style (dict): The complete style for plotting the dot.
    - applied_rule_label (str): The label to use for this dot.
    - legend_style_override (dict): The original style dict from the rule that provided the label.
    """
    final_style = {}
    base_style = {'color': '#808080', 'alpha': 0.3, 'size': config['dot_size'], 'edgecolor': None, 'linewidth': 0, 'label': 'Other', 'zorder': 1}
    final_style.update(base_style)
    applied_rule_label = base_style['label']

    deletion_num = None
    if isinstance(deletion_count, (int, float)) and not np.isnan(deletion_count):
        deletion_num = int(deletion_count)

    has_stop = '*' in str(sequence)
    sequence_str = str(sequence)
    # Assuming WT and correctly framed sequences end with 'AGVFQ'
    is_frameshift = not sequence_str.endswith('AGVFQ') if sequence_str else True


    # --- Step 1: Find the highest-priority BASE rule ---
    highest_priority_matched = float('inf')
    legend_style_override = None # <<< NEW: Will store the style dict of the rule that sets the label

    for i, rule in enumerate(BASE_STYLING_RULES):
        match = False
        conditions = rule['condition']
        style = rule['style']

        if not conditions: match = True # Default rule
        elif conditions.get('has_stop_codon') and has_stop: match = True
        elif conditions.get('is_frameshift') and is_frameshift: match = True
        elif 'value_range_x' in conditions and x_value is not None and not np.isnan(x_value):
            low, high = conditions['value_range_x']
            if (low is None or x_value >= low) and (high is None or x_value <= high): match = True
        elif 'value_range_y' in conditions and y_value is not None and not np.isnan(y_value):
            low, high = conditions['value_range_y']
            if (low is None or y_value >= low) and (high is None or y_value <= high): match = True
        elif 'deletion_range' in conditions and deletion_num is not None:
            low, high = conditions['deletion_range']
            if (low is None or deletion_num >= low) and (high is None or deletion_num <= high): match = True

        if match and i < highest_priority_matched:
            highest_priority_matched = i
            applied_rule_label = style.get('label', f"Rule {i}")
            if style.get('label'): # Store the base rule style if it has a label
                legend_style_override = style
            
            # Apply all styles from the base rule
            if 'color' in style and style['color'] is not None: final_style['color'] = style['color']
            if 'alpha' in style and style['alpha'] is not None: final_style['alpha'] = style['alpha']
            
            # Preserving original (buggy) size/zorder logic to minimize changes
            if 'size' in style and style['size'] is not None: final_style['size'] = style['size']
            if 'zorder' in style and style['zorder'] is not None: final_style['zorder'] = style['zorder']
            else: final_style['size'] = config['dot_size']

            if 'edgecolor' in style: final_style['edgecolor'] = style['edgecolor']
            if 'linewidth' in style: final_style['linewidth'] = style['linewidth'] if style['linewidth'] is not None else 0

    # --- Step 2: Apply all matching OVERLAY rules ---
    for i, rule in enumerate(OVERLAY_STYLING_RULES):
        match = False
        conditions = rule['condition']
        style = rule['style']
        
        # This loop's match logic only needs to check for overlay conditions
        if conditions.get('sequence_list') and sequence in conditions['sequence_list']: 
            match = True
        # (You could add other condition types here if needed)
        
        if match:
            # This rule matches, UPDATE final_style
            # This logic allows overlays to *not* specify a color, preserving the base color
            if style.get('label'): # Only update label if overlay provides one
                applied_rule_label = style.get('label')
                legend_style_override = style # <<< NEW: Store the rule's style dict
            
            if 'color' in style and style['color'] is not None: final_style['color'] = style['color']
            if 'alpha' in style and style['alpha'] is not None: final_style['alpha'] = style['alpha']
            
            # Original logic for overlay rules didn't have the buggy else
            if 'size' in style and style['size'] is not None: final_style['size'] = style['size']
            if 'zorder' in style and style['zorder'] is not None: final_style['zorder'] = style['zorder']
    
            if 'edgecolor' in style: final_style['edgecolor'] = style['edgecolor']
            if 'linewidth' in style: final_style['linewidth'] = style['linewidth'] if style['linewidth'] is not None else 0

    # --- Step 3: Finalize (same as before) ---
    if final_style.get('edgecolor') is None or final_style.get('linewidth') is None or final_style.get('linewidth') <= 0:
        final_style['edgecolor'] = 'none'
        final_style['linewidth'] = 0

    # Handle case where zorder/size might not have been set
    if 'size' not in final_style: final_style['size'] = config['dot_size']
    if 'zorder' not in final_style: final_style['zorder'] = 1 # Default zorder if not set by base


    return final_style, applied_rule_label, legend_style_override # <<< MODIFIED
# === EDIT END ===

# ==============================================================================
#                        PLOTTING FUNCTION
# ==============================================================================

def plot_replicate_comparison(file1_info, file2_info,
                              Frequency_df, log2_Frequency_df, fold_change_df, log2_fold_change_df,
                              deletion_df,
                              plot_config,
                              appearance_config,
                              fold_change_is_active, predefined_names_dict):
    """Generates a scatter plot comparing two replicates based on the provided plot_config (axis_type)."""
    start_time = time.time()

    f1_name = file1_info['full_name']
    f2_name = file2_info.get('full_name') if file2_info else None

    label1 = f"{file1_info['num']}{file1_info['rep_id']}{file1_info['label']}" if file1_info['rep_id'] else f"{file1_info['num']}_{file1_info['label']}"
    label2 = f"{file2_info['num']}{file2_info['rep_id']}{file2_info['label']}" if file2_info and file2_info['rep_id'] else (f"{file2_info['num']}_{file2_info['label']}" if file2_info else "N/A (Single Replicate)")

    current_axis_type = plot_config['axis_type']
    plot_title_type_suffix = current_axis_type

    effective_axis_type = current_axis_type
    if not fold_change_is_active:
        if current_axis_type == 'fold_change':
            print(f"  INFO: Fold Change not calculated. Plotting 'raw_Frequency' instead of 'fold_change' for {label1} vs {label2}.")
            effective_axis_type = 'raw_Frequency'
        elif current_axis_type == 'log2_fold_change':
            print(f"  INFO: Fold Change not calculated. Plotting 'log2_Frequency' instead of 'log2_fold_change' for {label1} vs {label2}.")
            effective_axis_type = 'log2_Frequency'

    plot_title_base = f"{label1} vs {label2}"
    print(f"Plotting: {label1} vs {label2} (Requested Axis: {current_axis_type}, Effective Axis: {effective_axis_type})")

    seqs1 = set(Frequency_df.index[Frequency_df[f1_name].notna()]) if f1_name in Frequency_df else set()
    seqs2 = set(Frequency_df.index[Frequency_df[f2_name].notna()]) if f2_name and f2_name in Frequency_df else set()
    all_sequences = list(seqs1.union(seqs2))

    if not all_sequences:
        print("  No sequence data found for either replicate. Skipping plot.")
        return

    data = {'sequence': all_sequences}
    data['deletion'] = [deletion_df.get(seq) for seq in all_sequences]

    current_axis_type_label_for_plot = ""

    if effective_axis_type == 'raw_Frequency':
        current_axis_type_label_for_plot = "Raw Frequency (%)"
        data['y'] = [Frequency_df.loc[seq, f1_name] if seq in seqs1 and f1_name in Frequency_df else 0 for seq in all_sequences]
        data['x'] = [Frequency_df.loc[seq, f2_name] if f2_name and seq in seqs2 and f2_name in Frequency_df else 0 for seq in all_sequences]
    elif effective_axis_type == 'log2_Frequency':
        current_axis_type_label_for_plot = ""
        data['y'] = [log2_Frequency_df.loc[seq, f1_name] if seq in log2_Frequency_df.index and f1_name in log2_Frequency_df else np.nan for seq in all_sequences]
        data['x'] = [log2_Frequency_df.loc[seq, f2_name] if f2_name and seq in log2_Frequency_df.index and f2_name in log2_Frequency_df else np.nan for seq in all_sequences]
    elif effective_axis_type == 'fold_change':
        current_axis_type_label_for_plot = "Fold Change"
        data['y'] = [fold_change_df.loc[seq, f1_name] if seq in fold_change_df.index and f1_name in fold_change_df else np.nan for seq in all_sequences]
        data['x'] = [fold_change_df.loc[seq, f2_name] if f2_name and seq in fold_change_df.index and f2_name in fold_change_df else np.nan for seq in all_sequences]
    elif effective_axis_type == 'log2_fold_change':
        current_axis_type_label_for_plot = ""
        data['y'] = [log2_fold_change_df.loc[seq, f1_name] if seq in log2_fold_change_df.index and f1_name in log2_fold_change_df else np.nan for seq in all_sequences]
        data['x'] = [log2_fold_change_df.loc[seq, f2_name] if f2_name and seq in log2_fold_change_df.index and f2_name in log2_fold_change_df else np.nan for seq in all_sequences]

    y_axis_label = f"{label1}"
    x_axis_label = f"{label2}"
    #plot_title = f"{plot_title_base}\n({current_axis_type_label_for_plot})" #This is for plotting the file titles too
    plot_title = f"{current_axis_type_label_for_plot}"

    plot_data = pd.DataFrame(data)
    plot_data.replace([np.inf, -np.inf], np.nan, inplace=True)

    if appearance_config['plot_top_n'] is not None and isinstance(appearance_config['plot_top_n'], int) and appearance_config['plot_top_n'] < len(plot_data):
        plot_data['avg_value_for_plot_top_n'] = plot_data[['x', 'y']].mean(axis=1).fillna(plot_data[['x', 'y']].min(axis=1).min() - 1)
        plot_data = plot_data.nlargest(appearance_config['plot_top_n'], 'avg_value_for_plot_top_n')
        print(f"  Plotting top {appearance_config['plot_top_n']} sequences based on average {effective_axis_type}.")

    plot_data_valid = plot_data.dropna(subset=['x', 'y'], how='all').copy()

    if plot_data_valid.empty:
        print("  No valid data points to plot after handling NaNs/Infinities. Skipping plot.")
        return

    # --- R-squared calculation ---
    r_square_value = None
    if appearance_config.get('SHOW_R_SQUARE_IN_LEGEND', False) and len(plot_data_valid) >= 2:
        # Drop any remaining NaNs just for the linear regression calculation
        clean_data_for_corr = plot_data_valid[['x', 'y']].dropna()
        if len(clean_data_for_corr) >= 2:
            slope, intercept, r_value, p_value, std_err = linregress(clean_data_for_corr['x'], clean_data_for_corr['y'])
            r_square_value = r_value**2

    # === EDIT START: Update apply call to get legend_override ===
    # Apply styles and sort by zorder for correct plotting order
    style_results_df = plot_data_valid.apply(
        lambda row: pd.Series(
            get_style_for_sequence(row['sequence'], row['x'], row['y'], row['deletion'], appearance_config),
            index=['style_dict', 'rule_label', 'legend_override']
        ), 
        axis=1
    )
    plot_data_valid = pd.concat([plot_data_valid, style_results_df], axis=1)
    
    # Extract zorder from the style_dict for sorting
    plot_data_valid['zorder'] = plot_data_valid['style_dict'].apply(lambda x: x.get('zorder', 1))
    plot_data_valid.sort_values(by='zorder', ascending=True, inplace=True)
    # === EDIT END ===


    top_n_details_for_legend = []
    plot_data_valid['is_top_n_highlight'] = False
    plot_data_valid['top_n_color'] = None

    if appearance_config.get('SHOW_TOP_SEQUENCES_IN_LEGEND', False) and appearance_config.get('LEGEND_TOP_N_SEQUENCES', 0) > 0:
        n_top = appearance_config['LEGEND_TOP_N_SEQUENCES']
        plot_data_valid['avg_xy_for_legend'] = plot_data_valid[['x', 'y']].mean(axis=1)
        plot_data_valid_sorted = plot_data_valid.assign(
            sort_val=plot_data_valid['avg_xy_for_legend'].fillna(-np.inf)
        ).sort_values(by='sort_val', ascending=False)
        top_n_df = plot_data_valid_sorted.head(n_top)
        
        available_colors = list(mcolors.TABLEAU_COLORS.values())
        custom_colors_more = ['#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6', '#bfef45']
        if len(available_colors) < n_top:
            available_colors.extend(custom_colors_more * ((n_top // len(custom_colors_more)) + 1))

        for i, idx in enumerate(top_n_df.index):
            seq = top_n_df.loc[idx, 'sequence']
            name = predefined_names_dict.get(seq, "")
            color = available_colors[i % len(available_colors)]
            avg_xy_val = top_n_df.loc[idx, 'avg_xy_for_legend']
            del_count = deletion_df.get(seq, 'N/A')
            plot_data_valid.loc[idx, 'is_top_n_highlight'] = True
            plot_data_valid.loc[idx, 'top_n_color'] = color
            top_n_details_for_legend.append({
                'seq': seq, 'name': name, 'color': color,
                'avg_xy_val': avg_xy_val, 'del_count': del_count
            })

    plot_styles = []
    unique_legend_entries = {}

    # === EDIT START: Update legend entry creation ===
    for idx, row in plot_data_valid.iterrows():
        # Get the results from the new columns
        style_dict = row['style_dict']
        rule_label = row['rule_label']
        legend_override = row['legend_override'] # This is the style dict from the rule

        if row['is_top_n_highlight']:
            style_dict['color'] = row['top_n_color']
            style_dict['size'] = appearance_config['dot_size'] * appearance_config.get('TOP_SEQUENCE_MARKER_SIZE_FACTOR', 1.2)
            style_dict['edgecolor'] = '#000000'
            style_dict['linewidth'] = 0.75
            style_dict['alpha'] = max(style_dict.get('alpha', 0.75), 0.75)
            
        elif rule_label not in unique_legend_entries:
            # Start with the style_dict from the data point (for base color)
            legend_style = {
                'color': style_dict.get('color', '#808080'), 
                'alpha': style_dict.get('alpha', 1.0),
                'edgecolor': style_dict.get('edgecolor', 'none'), 
                'linewidth': style_dict.get('linewidth', 0)
            }

            # Now, check if this rule has a legend_override
            # (it should, either from base or overlay rule)
            if legend_override:
                # Apply style properties from the rule dict
                if 'alpha' in legend_override:
                    legend_style['alpha'] = legend_override['alpha']
                if 'edgecolor' in legend_override:
                    legend_style['edgecolor'] = legend_override['edgecolor']
                if 'linewidth' in legend_override:
                    legend_style['linewidth'] = legend_override['linewidth']
                
                # THIS IS THE USER'S REQUEST:
                # Check for the special 'legend_color' property
                if 'legend_color' in legend_override:
                    legend_style['color'] = legend_override['legend_color'] # Override the color
            
            unique_legend_entries[rule_label] = legend_style
            
        plot_styles.append(style_dict) # Add the *plotting* style
    # === EDIT END ===

    colors = [s['color'] for s in plot_styles]
    alphas = [s['alpha'] for s in plot_styles]
    sizes = [s['size'] for s in plot_styles]
    edgecolors = [s['edgecolor'] for s in plot_styles]
    linewidths = [s['linewidth'] for s in plot_styles]

    plt.style.use('default')
    plt.rcParams['font.family'] = appearance_config['font_family']
    plt.rcParams['font.size'] = appearance_config['font_size']
    fig, ax = plt.subplots(figsize=appearance_config['fig_size'])

    ax.scatter(
        plot_data_valid['x'], plot_data_valid['y'],
        c=colors, alpha=alphas, s=sizes,
        edgecolors=edgecolors, linewidths=linewidths, marker='o'
    )

    ax.set_xlabel(x_axis_label)
    ax.set_ylabel(y_axis_label)
    ax.set_title(plot_title, fontsize=appearance_config['font_size'] + 2)
    ax.grid(True, linestyle='--', color='grey', alpha=0.25, linewidth=0.5)

    if appearance_config['x_limit']: ax.set_xlim(appearance_config['x_limit'])
    else:
        min_x_data, max_x_data = plot_data_valid['x'].min(), plot_data_valid['x'].max()
        if pd.notna(min_x_data) and pd.notna(max_x_data) :
            pad_x = (max_x_data - min_x_data) * 0.05 if (max_x_data - min_x_data) > 1e-9 else 0.5
            ax.set_xlim(min(0, min_x_data) - pad_x if min_x_data > -pad_x and min_x_data > 0 else min_x_data - pad_x,
                        max(0, max_x_data) + pad_x if max_x_data < pad_x and max_x_data < 0 else max_x_data + pad_x)

    if appearance_config['y_limit']: ax.set_ylim(appearance_config['y_limit'])
    else:
        min_y_data, max_y_data = plot_data_valid['y'].min(), plot_data_valid['y'].max()
        if pd.notna(min_y_data) and pd.notna(max_y_data):
            pad_y = (max_y_data - min_y_data) * 0.05 if (max_y_data - min_y_data) > 1e-9 else 0.5
            ax.set_ylim(min(0, min_y_data) - pad_y if min_y_data > -pad_y and min_y_data > 0 else min_y_data - pad_y,
                        max(0, max_y_data) + pad_y if max_y_data < pad_y and max_y_data < 0 else max_y_data + pad_y)

    legend_patches = []
    
    # --- Add R-squared to legend ---
    if r_square_value is not None:
        # The unicode for R-squared is U+00B2
        r_square_label = f"R\u00b2 = {r_square_value:.2f}"
        r_square_patch = mpatches.Patch(color='none', label=r_square_label)
        legend_patches.append(r_square_patch)
        # Add a spacer if there are other legend entries
        if unique_legend_entries or (appearance_config.get('SHOW_TOP_SEQUENCES_IN_LEGEND', False) and top_n_details_for_legend):
             legend_patches.append(mpatches.Patch(color='none', label='')) # Spacer
             
    # --- EDIT: Create a single combined list of all rule labels for ordering ---
    all_rule_labels_ordered = [r['style'].get('label') for r in BASE_STYLING_RULES if r['style'].get('label')]
    all_rule_labels_ordered.extend([r['style'].get('label') for r in OVERLAY_STYLING_RULES if r['style'].get('label')])

    sorted_legend_keys = sorted(unique_legend_entries.keys(), key=lambda k: all_rule_labels_ordered.index(k) if k in all_rule_labels_ordered else float('inf'))
    
    for label in sorted_legend_keys:
        style = unique_legend_entries[label]
        ec = style['edgecolor'] if style['edgecolor'] != 'none' else None
        lw = style['linewidth'] if ec is not None else 0
        patch = mpatches.Patch(facecolor=style['color'], alpha=style['alpha'], edgecolor=ec, linewidth=lw, label=label)
        legend_patches.append(patch)

    legend_title_main = ""

    if appearance_config.get('SHOW_TOP_SEQUENCES_IN_LEGEND', False) and top_n_details_for_legend:
        if legend_patches: legend_patches.append(mpatches.Patch(color='none', label='')) 

        title_patch_label = f"Top {len(top_n_details_for_legend)} by Avg. X-Y ({current_axis_type_label_for_plot}):"
        legend_patches.append(mpatches.Patch(color='none', label=title_patch_label))

        for detail in top_n_details_for_legend:
            seq_display_part = detail['seq']
            name_part = detail['name']
            del_part = f"Del: {detail['del_count']}"
            avg_val_num = detail['avg_xy_val']
            val_display_part = ""
            if pd.notna(avg_val_num):
                if "Raw Frequency (%)" in current_axis_type_label_for_plot: val_display_part = f"{avg_val_num:.2f}%"
                elif "Log2" in current_axis_type_label_for_plot: val_display_part = f"{avg_val_num:.2f}"
                elif "Fold Change" in current_axis_type_label_for_plot: val_display_part = f"{avg_val_num:.2f}x"
                else: val_display_part = f"{avg_val_num:.2f}"
            else: val_display_part = "N/A"

            full_label = seq_display_part
            if name_part: full_label += f" | {name_part}"
            full_label += f" | {del_part} | {val_display_part}"

            max_len = 60
            if len(full_label) > max_len:
                if name_part:
                    other_parts_len = len(f" | {name_part} | {del_part} | {val_display_part}")
                    available_seq_len = max_len - other_parts_len - 4
                    if available_seq_len > 5:
                        seq_display_part = detail['seq'][:available_seq_len] + "..."
                        full_label = f"{seq_display_part} | {name_part} | {del_part} | {val_display_part}"
                    else:
                        full_label = f"{name_part} | {del_part} | {val_display_part}"
                        if len(full_label) > max_len:
                           name_part_trunc = name_part[:max_len - len(f" | {del_part} | {val_display_part}") - 4] + "..."
                           full_label = f"{name_part_trunc} | {del_part} | {val_display_part}"
                else:
                    seq_display_part = detail['seq'][:max_len - len(f" | {del_part} | {val_display_part}") - 4] + "..."
                    full_label = f"{seq_display_part} | {del_part} | {val_display_part}"
            if len(full_label) > max_len:
                full_label = full_label[:max_len-3] + "..."

            patch = mpatches.Patch(color=detail['color'], label=full_label)
            legend_patches.append(patch)
        if legend_title_main == "Legend" and len(legend_patches) > len(top_n_details_for_legend)+2:
             legend_title_main = "Legend"
        else:
             legend_title_main = ""



    if legend_patches:
        # 1. Create the legend and store it in a variable 'leg'
        leg = ax.legend(
            handles=legend_patches, 
            title=legend_title_main, 
            loc='upper left', 
            bbox_to_anchor=(1.02, 1), 
            borderaxespad=0., 
            fontsize=appearance_config['font_size']-1
        )
        
        # 2. Get the legend's frame and set its properties
        frame = leg.get_frame()
        frame.set_edgecolor('grey')  # Or whatever color you want
        frame.set_linewidth(0.5)      # <--- EDIT THIS VALUE FOR LEGEND BORDER
        
        # --- This is your original code for the PLOT BORDER ---
        border_width = 0.5 # Set your desired width
        ax.spines['top'].set_linewidth(border_width)
        ax.spines['right'].set_linewidth(border_width)
        ax.spines['bottom'].set_linewidth(border_width)
        ax.spines['left'].set_linewidth(border_width)

        plt.tight_layout(rect=[0, 0, 0.80, 1])
    else:
        plt.tight_layout()

    output_base_filename = f"C_{label1}_{label2}_{plot_title_type_suffix}"
    output_base_filename = re.sub(r'[\\/*?:"<>|]', "_", output_base_filename)
    output_dir_path = Path(OUTPUT_DIRECTORY)
    output_dir_path.mkdir(parents=True, exist_ok=True)

    for fmt in appearance_config['save_formats']:
        output_path = output_dir_path / f"{output_base_filename}.{fmt}"
        try:
            plt.savefig(output_path, format=fmt, dpi=appearance_config['dpi'], bbox_inches='tight')
            print(f"  Plot saved as {fmt.upper()} at: {output_path}")
        except Exception as e:
            print(f"  Error saving {fmt.upper()}: {e}")

    plt.close(fig)
    end_time = time.time()
    print(f"  Plotting for {effective_axis_type} took {end_time - start_time:.2f} seconds.")

# ==============================================================================
#                            MAIN EXECUTION
# ==============================================================================

def main():
    """Main function to orchestrate the analysis and plotting."""
    overall_start_time = time.time()
    script_dir = Path(INPUT_DIRECTORY).resolve()
    excel_files = glob.glob(os.path.join(script_dir, "*.xlsx"))

    Path(OUTPUT_DIRECTORY).mkdir(parents=True, exist_ok=True)

    if not excel_files:
        print(f"Error: No Excel files (.xlsx) found in the directory: {script_dir}")
        return
    print(f"Found {len(excel_files)} Excel files in {script_dir}")

    file_data = {}
    all_sequences = set()
    raw_Frequency_dict = {}
    raw_reads_dict = {}
    deletion_dict = {}

    print("Processing files...")
    for file_path in excel_files:
        file_name_basename = os.path.basename(file_path)
        file_start_time = time.time()
        print(f"  Processing: {file_name_basename}...")
        num, rep_id, label = parse_filename(file_name_basename)
        if num is None:
            print(f"    Skipping file (could not parse name): {file_name_basename}")
            continue

        full_name_key = file_name_basename
        file_info = {'num': num, 'rep_id': rep_id, 'label': label, 'full_name': full_name_key, 'path': file_path}
        if num not in file_data: file_data[num] = []
        file_data[num].append(file_info)

        try:
            df = pd.read_excel(file_path)
            sequence_col_idx, reads_col_idx, Frequency_col_idx, deletion_col_idx = 1, 3, 4, 9

            if not all(idx < len(df.columns) for idx in [sequence_col_idx, reads_col_idx, Frequency_col_idx, deletion_col_idx]):
                print(f"    Skipping file: Missing required columns (B, D, E, or J). File has {len(df.columns)} columns.")
                if num in file_data:
                     file_data[num] = [f for f in file_data[num] if f['full_name'] != full_name_key]
                     if not file_data[num]: del file_data[num]
                continue

            sequence_col_name, reads_col_name, Frequency_col_name, deletion_col_name = df.columns[sequence_col_idx], df.columns[reads_col_idx], df.columns[Frequency_col_idx], df.columns[deletion_col_idx]

            for _, row in df.iterrows():
                sequence = str(row[sequence_col_name]).strip() if pd.notna(row[sequence_col_name]) else None
                if not sequence or sequence.lower() == 'nan': continue
                all_sequences.add(sequence)

                Frequency = pd.to_numeric(row[Frequency_col_name], errors='coerce')
                if pd.isna(Frequency): Frequency = 0.0

                reads = pd.to_numeric(row[reads_col_name], errors='coerce')
                if pd.isna(reads): reads = 0

                if sequence not in raw_Frequency_dict: raw_Frequency_dict[sequence] = {}
                raw_Frequency_dict[sequence][full_name_key] = Frequency
                if sequence not in raw_reads_dict: raw_reads_dict[sequence] = {}
                raw_reads_dict[sequence][full_name_key] = reads
                if sequence not in deletion_dict:
                    del_val = pd.to_numeric(row[deletion_col_name], errors='coerce')
                    deletion_dict[sequence] = int(del_val) if pd.notna(del_val) else np.nan
            print(f"    Processed {len(df)} rows in {time.time() - file_start_time:.2f} sec.")
        except Exception as e:
            print(f"    Error processing file {file_name_basename}: {e}")
            if num in file_data:
                 file_data[num] = [f for f in file_data[num] if f['full_name'] != full_name_key]
                 if not file_data[num]: del file_data[num]
            continue

    if not file_data or not any(file_data.values()):
        print("Error: No valid files could be processed successfully.")
        return
    print("File processing complete.")

    print("Creating DataFrames...")
    all_processed_files = sorted(list(set(fi['full_name'] for num_key in file_data for fi in file_data[num_key])))
    all_sequences = sorted(list(all_sequences))

    Frequency_df = pd.DataFrame.from_dict(raw_Frequency_dict, orient='index', columns=all_processed_files).reindex(index=all_sequences, columns=all_processed_files).fillna(0.0)
    reads_df = pd.DataFrame.from_dict(raw_reads_dict, orient='index', columns=all_processed_files).reindex(index=all_sequences, columns=all_processed_files).fillna(0)
    
    print("Calculating Log2 Raw Frequency...")
    log2_Frequency_df = Frequency_df.applymap(calculate_log2_value)

    # --- Fold Change Calculation using RPM ---
    fold_change_is_active = False
    fold_change_df = pd.DataFrame(index=all_sequences, columns=all_processed_files, dtype=float)
    log2_fold_change_df = pd.DataFrame(index=all_sequences, columns=all_processed_files, dtype=float)

    baseline_file_info = None
    baseline_file_name_key = None
    if 1 in file_data and file_data[1]:
        baseline_file_info = file_data[1][0]
        if len(file_data[1]) > 1:
            print(f"  WARNING: Multiple files for baseline condition 1. Using first: {baseline_file_info['full_name']}")
        baseline_file_name_key = baseline_file_info['full_name']
        fold_change_is_active = True
        print(f"Using baseline file: '{baseline_file_name_key}' for Fold Change calculation.")
    else:
        print("  INFO: No file for condition '1' found. Fold Change calculation will be skipped.")

    if fold_change_is_active:
        print("Calculating Fold Change and Log2 Fold Change using RPM method...")
        
        # 1. Calculate RPM for all samples
        total_reads_per_sample = reads_df.sum(axis=0)
        rpm_df = reads_df.apply(lambda col: ((col + 1) / total_reads_per_sample[col.name]) * 1_000_000 if total_reads_per_sample[col.name] > 0 else 0, axis=0)
        
        rpm_input_col = rpm_df[baseline_file_name_key]

        # 2. Iterate through selection files to calculate FC
        for sel_file_key in all_processed_files:
            if sel_file_key == baseline_file_name_key:
                fold_change_df[sel_file_key] = 1.0
                continue

            print(f"  Calculating FC for '{sel_file_key}'...")
            rpm_selection_col = rpm_df[sel_file_key]
            
            # Default FC calculation
            # Use .div and fill_value=np.inf for safe division where rpm_input_col could be 0
            fc_col = rpm_selection_col.div(rpm_input_col).replace(np.inf, np.nan)

            # Apply filter threshold logic if enabled
            if USE_FILTER_THRESHOLD:
                reads_selection_col = reads_df[sel_file_key]
                filter_rpm_value = None
                actual_threshold_used = None

                # Start from the defined threshold and search downwards
                for current_threshold in range(FILTER_READ_THRESHOLD, 0, -1):
                    sequences_at_threshold = reads_selection_col[reads_selection_col == current_threshold]
                    
                    if not sequences_at_threshold.empty:
                        rpms_at_threshold = rpm_df.loc[sequences_at_threshold.index, sel_file_key]
                        filter_rpm_value = rpms_at_threshold.max()
                        actual_threshold_used = current_threshold
                        break # Found the highest possible threshold, exit loop
                
                if filter_rpm_value is not None and filter_rpm_value > 0:
                    print(f"    Filter RPM for '{sel_file_key}' is: {filter_rpm_value:.2f} (determined from sequences with {actual_threshold_used} reads).")
                    
                    # Apply filter to sequences with reads LESS THAN the original threshold
                    low_reads_mask = reads_selection_col < FILTER_READ_THRESHOLD
                    
                    # Update FC only for the low-read sequences
                    fc_col[low_reads_mask] = rpm_selection_col[low_reads_mask] / filter_rpm_value
                else:
                    # This warning will now only appear if there are no sequences with reads between 1 and the threshold
                    print(f"    WARNING: No sequences found with reads between 1 and {FILTER_READ_THRESHOLD} for '{sel_file_key}'. Cannot apply filter logic. Using standard FC calculation for all sequences.")

            fold_change_df[sel_file_key] = fc_col

        # 3. Calculate Log2 Fold Change
        log2_fold_change_df = fold_change_df.applymap(calculate_log2_value)

    # --- WT Filtering ---
    if not INCLUDE_WT_IN_OUTPUT:
        print(f"\nExcluding WT sequence ('{WT_SEQUENCE}') from analysis and outputs.")
        if WT_SEQUENCE in all_sequences:
            # Use list comprehension to create a new list without WT
            all_sequences = [s for s in all_sequences if s != WT_SEQUENCE]
            Frequency_df.drop(WT_SEQUENCE, inplace=True, errors='ignore')
            log2_Frequency_df.drop(WT_SEQUENCE, inplace=True, errors='ignore')
            reads_df.drop(WT_SEQUENCE, inplace=True, errors='ignore')
            if fold_change_is_active:
                fold_change_df.drop(WT_SEQUENCE, inplace=True, errors='ignore')
                log2_fold_change_df.drop(WT_SEQUENCE, inplace=True, errors='ignore')
            if WT_SEQUENCE in deletion_dict:
                del deletion_dict[WT_SEQUENCE]
        else:
            print(f"  Warning: WT sequence '{WT_SEQUENCE}' not found in the dataset.")

    # --- NEW: Deletion Count Filtering ---
    if FILTER_BY_DELETION_COUNT:
        print(f"\nFiltering sequences with deletions <= {DELETION_COUNT_THRESHOLD}...")
        
        # Identify sequences to remove based on the deletion count
        sequences_to_remove = {
            seq for seq, del_count in deletion_dict.items()
            if pd.notna(del_count) and del_count <= DELETION_COUNT_THRESHOLD
        }
        print(f"  Found {len(sequences_to_remove)} sequences that meet the deletion criteria.")

        # Ensure WT is not in the removal list
        if WT_SEQUENCE in sequences_to_remove:
            sequences_to_remove.remove(WT_SEQUENCE)
            print(f"  Exempting WT sequence ('{WT_SEQUENCE}') from this filter.")

        if sequences_to_remove:
            print(f"  Removing {len(sequences_to_remove)} sequences from all datasets.")
            sequences_to_remove_list = list(sequences_to_remove)
            
            # Update all relevant data structures
            all_sequences = [s for s in all_sequences if s not in sequences_to_remove_list]
            Frequency_df.drop(index=sequences_to_remove_list, inplace=True, errors='ignore')
            log2_Frequency_df.drop(index=sequences_to_remove_list, inplace=True, errors='ignore')
            reads_df.drop(index=sequences_to_remove_list, inplace=True, errors='ignore')
            if fold_change_is_active:
                fold_change_df.drop(index=sequences_to_remove_list, inplace=True, errors='ignore')
                log2_fold_change_df.drop(index=sequences_to_remove_list, inplace=True, errors='ignore')
            
            for seq in sequences_to_remove_list:
                if seq in deletion_dict:
                    del deletion_dict[seq]
        else:
            print("  No sequences to remove after applying the filter.")
            
    plot_types_to_generate = []
    if PLOT_RAW_Frequency: plot_types_to_generate.append('raw_Frequency')
    if PLOT_LOG2_Frequency: plot_types_to_generate.append('log2_Frequency')
    if PLOT_FOLD_CHANGE: plot_types_to_generate.append('fold_change')
    if PLOT_LOG2_FOLD_CHANGE: plot_types_to_generate.append('log2_fold_change')

    if plot_types_to_generate:
        print(f"\nGenerating plots for types: {', '.join(plot_types_to_generate)}...")
        conditions_to_plot_nums = sorted(file_data.keys())

        for plot_axis_type in plot_types_to_generate:
            current_plot_config = {'axis_type': plot_axis_type} 

            print(f"\n--- Generating plots for axis type: {plot_axis_type} ---")
            for num in conditions_to_plot_nums:
                replicate_files_info = file_data[num]
                
                if fold_change_is_active and num == 1 and len(file_data) > 1:
                     if len(replicate_files_info) <= 1:
                         print(f"  Skipping plots for Condition 1 as it's the baseline and has no internal replicates.")
                         continue
                     else:
                         print(f"  Condition 1 is baseline, plotting comparisons within its replicates.")

                if len(replicate_files_info) == 1:
                    plot_replicate_comparison(
                        file1_info=replicate_files_info[0], file2_info=None,
                        Frequency_df=Frequency_df, log2_Frequency_df=log2_Frequency_df,
                        fold_change_df=fold_change_df, log2_fold_change_df=log2_fold_change_df,
                        deletion_df=deletion_dict, plot_config=current_plot_config,
                        appearance_config=PLOT_APPEARANCE_CONFIG, fold_change_is_active=fold_change_is_active,
                        predefined_names_dict=predefined_names)
                elif len(replicate_files_info) > 1:
                    for file1_info, file2_info in itertools.combinations(replicate_files_info, 2):
                        plot_replicate_comparison(
                            file1_info=file1_info, file2_info=file2_info,
                            Frequency_df=Frequency_df, log2_Frequency_df=log2_Frequency_df,
                            fold_change_df=fold_change_df, log2_fold_change_df=log2_fold_change_df,
                            deletion_df=deletion_dict, plot_config=current_plot_config,
                            appearance_config=PLOT_APPEARANCE_CONFIG, fold_change_is_active=fold_change_is_active,
                            predefined_names_dict=predefined_names)
    else:
        print("\nNo plot types selected to generate.")

    if CREATE_EXCEL_OUTPUT_SUMMARY:
        print("\nCreating Excel Output Summaries (per comparison)...")
        conditions_for_excel_iteration = sorted(file_data.keys())

        for num_condition_key in conditions_for_excel_iteration:
            replicate_files_info_list = file_data[num_condition_key]

            def get_condensed_filename_part(file_info_dict):
                num_str, rep_id_str, label_str = str(file_info_dict['num']), file_info_dict['rep_id'] or "", file_info_dict['label']
                return f"{num_str}{label_str}{rep_id_str}"

            if not replicate_files_info_list: continue

            if len(replicate_files_info_list) == 1:
                file1_info = replicate_files_info_list[0]
                f1_key = file1_info['full_name']
                part1_name = get_condensed_filename_part(file1_info)
                excel_filename_base = f"C_{part1_name}_vs_Single"
                excel_filename_base = re.sub(r'[\\/*?:"<>|]', "_", excel_filename_base)
                excel_full_filename = f"{excel_filename_base}.xlsx"
                
                print(f"    Preparing summary: {excel_full_filename}")
                if f1_key not in Frequency_df.columns: continue
                
                files_to_average_keys = [f1_key]

                current_summary_data = []
                for seq in all_sequences:
                    avg_fc, avg_log2_fc = np.nan, np.nan
                    if fold_change_is_active:
                        avg_fc = fold_change_df.loc[seq, files_to_average_keys].mean()
                        avg_log2_fc = log2_fold_change_df.loc[seq, files_to_average_keys].replace([-np.inf, np.inf], np.nan).mean()
                    current_summary_data.append({
                        'Sequence': seq, 'Name': predefined_names.get(seq, ""), 'Deletions': deletion_dict.get(seq, np.nan),
                        'Average Reads': reads_df.loc[seq, files_to_average_keys].mean(),
                        'Average %Frequency': Frequency_df.loc[seq, files_to_average_keys].mean(),
                        'Average log2 %Frequency': log2_Frequency_df.loc[seq, files_to_average_keys].replace([-np.inf, np.inf], np.nan).mean(),
                        'Average Fold Change': avg_fc, 'Average Log2 Fold Change': avg_log2_fc})
                
                if not current_summary_data: continue
                summary_df = pd.DataFrame(current_summary_data)
                sort_cols = ['Average Log2 Fold Change', 'Average Fold Change', 'Average log2 %Frequency']
                for col in sort_cols:
                    if col in summary_df.columns and not summary_df[col].isnull().all():
                        summary_df = summary_df.sort_values(by=col, ascending=False).reset_index(drop=True)
                        break
                summary_df.insert(0, 'Rank', summary_df.index + 1)
                summary_file_path = Path(OUTPUT_DIRECTORY) / excel_full_filename
                summary_df.to_excel(summary_file_path, index=False)
                print(f"      Summary saved to: {summary_file_path}")
            
            elif len(replicate_files_info_list) > 1:
                for file1_info, file2_info in itertools.combinations(replicate_files_info_list, 2):
                    f1_key, f2_key = file1_info['full_name'], file2_info['full_name']
                    part1_name, part2_name = get_condensed_filename_part(file1_info), get_condensed_filename_part(file2_info)
                    excel_filename_base = f"C_{part1_name}_{part2_name}"
                    excel_filename_base = re.sub(r'[\\/*?:"<>|]', "_", excel_filename_base)
                    excel_full_filename = f"{excel_filename_base}.xlsx"

                    print(f"    Preparing summary: {excel_full_filename}")
                    if f1_key not in Frequency_df.columns or f2_key not in Frequency_df.columns: continue
                    
                    files_to_average_keys = [f1_key, f2_key]
                    
                    current_summary_data = []
                    for seq in all_sequences:
                        avg_fc, avg_log2_fc = np.nan, np.nan
                        if fold_change_is_active:
                            avg_fc = fold_change_df.loc[seq, files_to_average_keys].mean()
                            avg_log2_fc = log2_fold_change_df.loc[seq, files_to_average_keys].replace([-np.inf, np.inf], np.nan).mean()
                        current_summary_data.append({
                            'Sequence': seq, 'Name': predefined_names.get(seq, ""), 'Deletions': deletion_dict.get(seq, np.nan),
                            'Average Reads': reads_df.loc[seq, files_to_average_keys].mean(),
                            'Average %Frequency': Frequency_df.loc[seq, files_to_average_keys].mean(),
                            'Average log2 %Frequency': log2_Frequency_df.loc[seq, files_to_average_keys].replace([-np.inf, np.inf], np.nan).mean(),
                            'Average Fold Change': avg_fc, 'Average Log2 Fold Change': avg_log2_fc})

                    if not current_summary_data: continue
                    summary_df = pd.DataFrame(current_summary_data)
                    sort_cols = ['Average Log2 Fold Change', 'Average Fold Change', 'Average log2 %Frequency']
                    for col in sort_cols:
                        if col in summary_df.columns and not summary_df[col].isnull().all():
                            summary_df = summary_df.sort_values(by=col, ascending=False).reset_index(drop=True)
                            break
                    summary_df.insert(0, 'Rank', summary_df.index + 1)
                    summary_file_path = Path(OUTPUT_DIRECTORY) / excel_full_filename
                    summary_df.to_excel(summary_file_path, index=False)
                    print(f"      Summary saved to: {summary_file_path}")
    else:
        print("\nExcel output summary creation is disabled.")

    overall_end_time = time.time()
    print(f"\nScript finished. Total execution time: {overall_end_time - overall_start_time:.2f} seconds.")
    warnings.filterwarnings("default")

if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=RuntimeWarning, message="Mean of empty slice")
    warnings.filterwarnings("ignore", category=UserWarning, message="This figure includes Axes that are not compatible with tight_layout")
    main()

