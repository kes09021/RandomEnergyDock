# Random Energy Dock

This repository contains various scripts for protein and ligand docking.

## Requirements
The following Python packages are required to run these scripts:
- pymol
- biopython
- rdkit
- openbabel
- scikit-learn
- numpy

Installation:
pip install -r requirements.txt

## Script Descriptions

### cif_to_pdbqt.py
Converts protein CIF files to PDBQT files.
Usage:
python cif_to_pdbqt.py <CIF file path>

### pdb_to_pdbqt.py
Converts protein PDB files to PDBQT files.
Usage:
python pdb_to_pdbqt.py <PDB file path>

### smiles_to_pdbqt.py
Converts SMILES format ligands to PDBQT files.
Usage:
python smiles_to_pdbqt.py <SMILES>

### random_docking.py
Performs random docking using ligand and protein PDB and PDBQT files.
Usage:
python random_docking.py <protein_pdbqt> <ligand_pdbqt> <protein_pdb>

### visualize_energy_pymol.py
Calculates binding free energy from docking results and visualizes it using PyMOL.
Usage:
python visualize_energy_pymol.py <protein PDB file path> <PDBQT file path/docking_output_*.pdbqt>

## Visualization
The following image is an example of a protein-ligand docking visualization created by this tool. The red dots represent regions where the binding free energy is stable, indicating potential binding sites.


![Protein-Ligand Docking Visualization](path/to/image.png)