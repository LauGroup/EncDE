# Directed evolution of multimeric proteins is enabled by dual-compensatory gene duplication

## Overview
This repository contains the custom Python scripts and bioinformatics workflow used to analyze Directed Evolution (DE) experiments described in our manuscript: **"Directed evolution of multimeric proteins is enabled by dual-compensatory gene duplication"**.

The workflow processes Next-Generation Sequencing (NGS) data to track the evolution of protein variants, focusing on the enrichment of unique sequences.

## Repository Organization
* **`Scripts/`**: Python scripts for the analysis pipeline including Fastq processing, alignment, and visualization.
* **`Gen_1_LR/` to `Gen_4_Mixed/`**: Experimental datasets containing processed data, alignments, and output plots for different generations of the directed evolution experiments.
* **`Logos/`**: Sequence logo visualizations generated from the analysis.

## Dependencies
The scripts are designed to run in a standard Python environment. The alignment step requires Clustal Omega.

* **Python 3.x**
* **Required Python Packages:**
    * `pandas`
    * `biopython`
    * `logomaker`
    * `matplotlib`
* **External Tools:**
    * **Clustal Omega**: Binaries for macOS and Windows are included in `Scripts/clustalo mac/` and `Scripts/clustalo windows/`.

## Workflow Description
The analysis pipeline consists of the following main stages:

1.  **`1_Fastq_to_Excel.py`**: Converts raw FASTQ sequencing reads into frequency count tables.
2.  **`2_Excel_to_Aln_ClustalO...py`**: Performs multiple sequence alignment (MSA) using Clustal Omega to identify mutations and deletions. There are two versions of this script for MAC or WINDOWS users.
3.  **`3_Aln_to_CountDels.py`**: Analyzes alignments to quantify the frequency of specific deletion patterns.
4.  **`4_Plot_replicates.py`**: Calculates frequency across replicates and fold-change over input library across replicates. Generates scatter plots and summaries in excel format.
5.  **`5_Aln_to_Logomaker.py`**: Generates sequence logos of frequency count tables to show distribution of sequences.

---
## How to analyse the data
This workflow will take you from the raw fastq files to final plots for analysis. 

All fastq files are found in NCBI Sequence Read Archive (SRA) database under the BioProject: [PRJNA1416590](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1416590)

Place fastq files in folders and define the directory in each script to process them. 

The analysed data is already available here in each folder **`Gen_1_LR/` to `Gen_4_Mixed/`**

## Demo: Processing the `Gen_2_MO` Dataset as a demo

The **Generation 2 (MO)** dataset is provided as a demo. 

The folder `Demo/Fastq/` contains raw fastq files to serve as inputs for processing.

The folder `Demo/` serves as the output for processed data. 

Open the repo folder on your interpreter like VScode

### Step 1: Fastq to Excel
1.  Open the script: `1_Fastq_to_Excel.py`.
2.  Set **Input:** `Demo/Fastq/` This folder contains the .fastq files to process. 
3.  Set **Output:** `Demo/` The script will generate frequency table excel files in the this directory.

Result: The script will parse through the fastq files according to the parameters set in the script, and will cluster reads for each unique sequence and create frequency tables. 

### Step 2: Excel to Alignment
1.  Run the script appropriate for your OS: `2_Excel_to_Aln_ClustalO_MAC.py` or `2_Excel_to_Aln_ClustalO_WINDOWS.py`.
2.  Set **Input:** `Demo/` This folder contains the excel files generated from the previous step.
3.  Set **Output:** `Demo/` The script will generate alignment files in the same directory by overwriting the excel files.

### Step 3: Analyze Deletions
1.  Run `python 3_Aln_to_CountDels.py`.
2.  Set **Input:** `Demo/` This folder contains the excel files with alignments generated from the previous step.
3.  Set **Output:** `Demo/` The script will overwrite the excel files from the previous step with a new column at the end with the number of deleted residues. 

### Step 4: Calculate fold change over input library across replicates and visualize results
1.  Open `python 4_Plot_replicates.py`.
2.  Set **Input:** `Demo/` 

The script decides file types according to the prefix numbers. 

The input libraries are files with "1" as prefix. (e.g., `1MO_Input_rep1.xlsx` and `1MO_Input_rep2.xlsx`).

The files to analyse have prefixes other than '1' and same prefix numbers indicate replicates. (e.g., `2MO_ExpComp_rep1.xlsx` and `2MO_ExpComp_rep2.xlsx`).

`ExpComp` means `pExplorer` and `pCompensator` was used during selection. 

`Exp` means only `pExplorer` was used during selection. 

3.  Set **Output:** `Demo/` The script will calculate frequency and fold change across replicates. Scatter plots showing the correlation between replicates and summaries in excel format will be generated. 

---

## Citing this workflow
Please cite our manuscript:

Siddiquee, R. et al. (2026) [*Directed evolution of multimeric proteins is enabled by dual-compensatory gene duplication.*](https://doi.org/10.64898/2026.01.12.698938) [Preprint]. doi:10.64898/2026.01.12.698938.

## References

Sievers F et al. (2011) [Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega](https://doi.org/10.1038/msb.2011.75) Mol Syst Biol. 7:539.

Tareen A & Kinney JB. (2019) [Logomaker: beautiful sequence logos in Python](https://doi.org/10.1093/bioinformatics/btz921) Bioinformatics. 36(7):2272–2274.