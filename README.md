# MS/MS Viewer

## Overview
**MS/MS Viewer** is a Python tool for visualizing peptide fragmentation spectra from mzXML files. Given an **mzXML file**, a **scan number**, and a **peptide sequence**, the program extracts the corresponding MS/MS spectrum, computes **b-ion** and **y-ion** m/z values, and overlays them onto the spectrum. It annotates matching peaks, helping users assess peptide-spectrum matches for peptide identification.

## Features
- Parses **mzXML** files to extract MS/MS spectra  
- Computes **b-ion** and **y-ion** fragment m/z values  
- Matches theoretical ion peaks to experimental data  
- Generates a **clear, labeled plot** for visualization  
- Helps determine if a peptide matches a given spectrum  

## Installation
Ensure you have Python installed, then install dependencies:
```sh
pip install matplotlib
```

## Usage
Run the script from the command line:
```sh
python script.py <mzxml_gzip_file> <scan_number> <peptide_sequence>
```
Example:
```sh
python script.py sample.mzXML.gz 1200 PEPTIDESEQ
```

## Output
- A **plot** displaying the MS/MS spectrum  
- Annotated peaks for **b-ion** and **y-ion** matches  
- A visual representation to assess peptide-spectrum matching  
