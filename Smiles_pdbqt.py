# smiles_to_pdbqt.py

import os
import sys

def smiles_to_pdb(smiles, pdb_file):
    os.system(f'obabel -:"{smiles}" -O {pdb_file} --gen3d')
    print(f'SMILES 문자열이 {pdb_file} 파일로 변환되었습니다.')

def pdb_to_pdbqt_ligand(pdb_file, pdbqt_file):
    os.system(f'obabel {pdb_file} -O {pdbqt_file} --partialcharge gasteiger')
    print(f'{pdb_file} 파일이 {pdbqt_file} 파일로 변환되었습니다.')

def main(smiles):
    ligand_pdb_file = "ligand.pdb"
    ligand_pdbqt_file = "ligand.pdbqt"

    smiles_to_pdb(smiles, ligand_pdb_file)
    pdb_to_pdbqt_ligand(ligand_pdb_file, ligand_pdbqt_file)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("사용 방법: python smiles_to_pdbqt.py <SMILES>")
        sys.exit(1)

    smiles = sys.argv[1]
    main(smiles)
