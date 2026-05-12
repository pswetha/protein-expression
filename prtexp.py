# ============================================================
# Computational Analysis of Human IL-2 (Target-X)
# Sequence + Structural Aggregation Risk + PyMOL Visualization
# ============================================================

import os
import csv
import subprocess
from Bio.Seq import Seq
from Bio.SeqUtils import ProtParam
from Bio.PDB import PDBParser, DSSP, NeighborSearch
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import numpy as np

# ------------------------------------------------------------
# Input sequence
# ------------------------------------------------------------
sequence = Seq(
    "MYRMQLLSCIALSLALVTNSAPTSSSTKKTQLQLEHLLLDLQMVILNGINNYKNPKLTRML"
    "TFKFYMPKKATELKHLQCLEEELKPLEEVLNLAQSKNFHLRPRDLISNINVIVLELKGSE"
    "TTFMCEYADEKTATIVEFLNRWITFCQSIISTLT"
)

PDB_FILE = "/home/vm-1/1M47.pdb"
CSV_OUT = "hydrophobic_surface_patches.csv"

HYDROPHOBIC = {"A", "V", "I", "L", "M", "F", "W", "Y"}
RSA_THRESHOLD = 0.25     # >25% exposed
DISTANCE_CUTOFF = 6.0    # Å


PYMOL_SCRIPT = "auto_color_aggregation.py"
PNG_FILE = "IL2_aggregation_patches.png"

# ------------------------------------------------------------
# 1. Physicochemical Profiling
# ------------------------------------------------------------
analyzer = ProtParam.ProteinAnalysis(str(sequence))

print("\n=== Physicochemical Properties ===")
print(f"Molecular Weight (Da): {analyzer.molecular_weight():.2f}")
print(f"Isoelectric Point (pI): {analyzer.isoelectric_point():.2f}")
print(f"GRAVY Score: {analyzer.gravy():.3f}")
print(f"Instability Index: {analyzer.instability_index():.2f}")

# ------------------------------------------------------------
# 2. Structural Evaluation
# ------------------------------------------------------------
# -----------------------------
# LOAD STRUCTURE
# -----------------------------
parser = PDBParser(QUIET=True)
structure = parser.get_structure("MODEL", PDB_FILE)
model = structure[0]

# -----------------------------
# DSSP
# -----------------------------
dssp = DSSP(model, PDB_FILE)

# -----------------------------
# Collect exposed hydrophobic residues
# -----------------------------
candidates = []

for (chain_id, res_id), d in dssp.property_dict.items():
    aa = d[1]
    rsa = d[3]

    if aa in HYDROPHOBIC and rsa >= RSA_THRESHOLD:
        chain = model[chain_id]
        residue = chain[res_id]

        if "CA" in residue:
            ca = residue["CA"]
            candidates.append({
            "chain": chain_id,
            "resname": residue.resname,
            "resseq": res_id[1],
            "aa": aa,
            "rsa": rsa,
            "atom": ca      # <-- store Atom object
            })


# -----------------------------
# Neighbor clustering
# -----------------------------
ns = NeighborSearch([c["atom"] for c in candidates])

patches = []

for r in candidates:
    neighbors = ns.search(
        r["atom"].get_coord(),
        DISTANCE_CUTOFF
    )

    if len(neighbors) >= 2:
        patches.append({
            "chain": r["chain"],
            "resname": r["resname"],
            "resseq": r["resseq"],
            "aa": r["aa"],
            "rsa": round(r["rsa"], 3),
            "neighbor_count": len(neighbors)
        })

df = pd.DataFrame(patches)

if df.empty:
    print("[!] No aggregation patches detected")
    print("[!] Try DISTANCE_CUTOFF = 7–8 Å or neighbor ≥ 2")
else:
    df.sort_values(by="neighbor_count", ascending=False, inplace=True)
    df.to_csv(CSV_OUT, index=False)
    print(f"[✔] Found {len(df)} aggregation-prone surface residues")

# ------------------------------------------------------------
# 5. Auto-generate PyMOL coloring script
# ------------------------------------------------------------
parser = PDBParser(QUIET=True)
structure = parser.get_structure("MODEL", PDB_FILE)
model = structure[0]

# Load aggregation residues
df = pd.read_csv(CSV_OUT)
agg_residues = set(df["resseq"].tolist())

# Collect coordinates
xs, ys, zs = [], [], []
colors = []

for chain in model:
    for res in chain:
        if "CA" in res:
            x, y, z = res["CA"].coord
            xs.append(x)
            ys.append(y)
            zs.append(z)

            # Color aggregation-prone residues
            if res.id[1] in agg_residues:
                colors.append("red")
            else:
                colors.append("lightgrey")

# Plot
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection="3d")

ax.scatter(xs, ys, zs, c=colors, s=25, alpha=0.9)

ax.set_title("IL-2 Aggregation-Prone Surface Residues", fontsize=12)
ax.set_axis_off()

plt.tight_layout()
plt.savefig(PNG_FILE, dpi=300)
plt.close()

print(f"[✔] Aggregation visualization saved as {PNG_FILE}")

print("\n=== Purification Recommendation ===")
print("Suggested buffer: 20 mM HEPES, 150 mM NaCl, pH 7.5")
print("Chromatography: Ni-NTA (IMAC) followed by Size Exclusion Chromatography")

