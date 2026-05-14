
# Protein Expression & Structural Aggregation Risk Analysis

Computational analysis pipeline for **Human IL-2 (Target-X)** protein expression profiling, structural aggregation prediction, and visualization using **Python + BioPython + PyMOL-compatible outputs**.

This project performs:

- **Physicochemical profiling**
- **Protein stability estimation**
- **Hydrophobic surface patch detection**
- **Aggregation-prone residue identification**
- **3D structural visualization**
- **Purification strategy recommendation**

---

## Project Overview

Protein aggregation is one of the major challenges in recombinant protein expression and purification.

This workflow identifies:

- Exposed hydrophobic residues
- Surface aggregation hotspots
- Structural instability regions
- Aggregation-prone clusters

using structural and sequence-level computational analysis.

---

## Features

### 1. Physicochemical Analysis
Calculates:

- Molecular Weight
- Isoelectric Point (pI)
- GRAVY Score
- Instability Index

---

### 2. Structural Analysis

Uses DSSP + PDB structure to detect:

- Solvent accessibility
- Hydrophobic exposure
- Surface aggregation patches

Parameters:

```python
RSA_THRESHOLD = 0.25
DISTANCE_CUTOFF = 6.0 Å
```

---

### 3. Aggregation Patch Detection

Identifies clustered exposed hydrophobic residues:

- Alanine (A)
- Valine (V)
- Isoleucine (I)
- Leucine (L)
- Methionine (M)
- Phenylalanine (F)
- Tryptophan (W)
- Tyrosine (Y)

Outputs:

```csv
hydrophobic_surface_patches.csv
```

---

### 4. 3D Visualization

Generates structural plot highlighting:

- **Red:** Aggregation-prone residues
- **Grey:** Non-risk residues

Output:

```png
IL2_aggregation_patches.png
```

---

### 5. Purification Recommendation

Provides suggested purification conditions:

- **Buffer:** 20 mM HEPES
- **Salt:** 150 mM NaCl
- **pH:** 7.5

Recommended chromatography:

- Ni-NTA (IMAC)
- Size Exclusion Chromatography (SEC)

---

# Repository Structure

```bash
.
├── prtexp.py
├── hydrophobic_surface_patches.csv
├── IL2_aggregation_patches.png
├── auto_color_aggregation.py
└── README.md
```

---

# Installation

## Clone Repository

```bash
git clone https://github.com/pswetha/protein-expression-analysis.git
cd protein-expression-analysis
```

---

## Install Dependencies

```bash
pip install biopython pandas numpy matplotlib
```

Install DSSP:

Ubuntu/Linux:

```bash
sudo apt-get install dssp
```

Install PyMOL (optional):

```bash
sudo apt-get install pymol
```

---

# Input Requirements

Edit input sequence in:

```python
sequence = Seq("YOUR_PROTEIN_SEQUENCE")
```

Set your structure path:

```python
PDB_FILE = "/path/to/your_structure.pdb"
```

---

# Usage

Run:

```bash
python prtexp.py
```

---

# Example Output

## Console Output

```bash
=== Physicochemical Properties ===
Molecular Weight: XXXX Da
Isoelectric Point: X.XX
GRAVY Score: X.XXX
Instability Index: XX.XX

[✔] Found aggregation-prone residues
[✔] Aggregation visualization saved
```

---

# Applications

Useful for:

- Recombinant protein production
- Therapeutic protein engineering
- Biopharmaceutical optimization
- Aggregation hotspot prediction
- Protein formulation design

---

# Technologies Used

- Python
- BioPython
- DSSP
- NumPy
- Pandas
- Matplotlib
- PyMOL

---

# Scientific Relevance

This workflow supports:

- Protein developability assessment
- Structural bioinformatics analysis
- Early-stage biologics screening
- Rational protein engineering

---

# Future Improvements

Planned additions:

- Machine learning aggregation prediction
- Rosetta integration
- Molecular dynamics simulation
- Interactive PyMOL automation
- Multi-protein batch analysis

---



**Pulakuntla Swetha**  
Protein Expression and Aggregation Risk Analysis Pipeline (2026)
