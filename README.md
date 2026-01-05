# biotech-gui-tools
Python-based GUI application for laboratory calculations (DNA analysis, Tm calculation, and molarity)
# BioTech Lab Assistant (GUI)

# BioTech OS Pro

**A comprehensive Bioinformatics Analysis Environment built with Streamlit and Python.**

![Python](https://img.shields.io/badge/Python-3.9+-blue?style=for-the-badge&logo=python)
![Streamlit](https://img.shields.io/badge/Streamlit-App-FF4B4B?style=for-the-badge&logo=streamlit)

##  Core Modules

* **Sequence Analysis:** DNA/RNA to Protein translation with ORF finder.
* **Restriction Mapping:** Virtual digestion with common enzymes (EcoRI, HindIII, etc.).
* **Protein Biophysics:** Molecular weight, Hydropathy profile, and Isoelectric Point (pI) calculation.
* **Alignment Engine:** Pairwise sequence alignment using the **Needleman-Wunsch algorithm**.
* **Genome Stats:** GC-skew analysis and nucleotide distribution tracking.

A desktop application designed to streamline routine laboratory calculations for biotechnologists. This tool provides a user-friendly interface for DNA analysis and PCR primer preparation.

## Current Features
* **Tm Calculation:** Estimates the melting temperature (Tm) of oligonucleotide primers using the **Wallace Rule** (Salt-adjusted formula for sequences <14-20 bp).
* **GC-Content Analysis:** Calculates the percentage of Guanine and Cytosine to predict primer stability and secondary structures.
* **Input Validation:** Built-in error handling to ensure only valid genomic data (A, T, G, C) is processed.

## Scientific Background
The melting temperature is calculated based on the following thermodynamic approximation:
$Tm = 2(A+T) + 4(G+C)$

This is a fundamental calculation for determining the optimal annealing temperature in PCR protocols.



## Technology Stack
* **Language:** Python 3.x
* **Library:** Tkinter (Standard GUI Library)

## Future Updates
* [ ] Molarity and dilution calculator.
* [ ] Reverse complement sequence generator.
* [ ] Protein molecular weight estimator.

## How to Run
1. Ensure you have Python installed.
2. Download `dna_calculator.py`.
3. Run the script: `python dna_calculator.py`
