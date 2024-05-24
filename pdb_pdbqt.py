# pdb_to_pdbqt.py

import os
import sys
from rdkit import Chem
from rdkit.Chem import AllChem

def pdb_to_pdbqt_protein(pdb_file, pdbqt_file):
    mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)

    if mol is None:
        print(f'{pdb_file} 파일을 로드하는 데 실패했습니다.')
        return

    if mol.GetNumConformers() == 0:
        AllChem.EmbedMolecule(mol)
        AllChem.UFFOptimizeMolecule(mol)

    AllChem.ComputeGasteigerCharges(mol)

    with open(pdbqt_file, 'w') as f:
        for atom in mol.GetAtoms():
            pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
            charge = atom.GetDoubleProp("_GasteigerCharge")
            f.write(f"ATOM  {atom.GetIdx()+1:>5}  {atom.GetSymbol():<2}  PRO A   1    {pos.x:>8.3f}{pos.y:>8.3f}{pos.z:>8.3f}  1.00  {charge:>6.2f}           {atom.GetSymbol():>2}\n")
        f.write("END\n")
    print(f'{pdb_file} 파일이 {pdbqt_file} 파일로 변환되었습니다.')

def main(pdb_file):
    protein_pdbqt_file = "protein.pdbqt"
    pdb_to_pdbqt_protein(pdb_file, protein_pdbqt_file)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("사용 방법: python pdb_to_pdbqt.py <PDB 파일 경로>")
        sys.exit(1)

    pdb_file = sys.argv[1]
    main(pdb_file)
