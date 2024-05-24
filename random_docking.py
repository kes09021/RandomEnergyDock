# random_docking.py

import os
import sys
import numpy as np
from Bio import PDB
from sklearn.cluster import KMeans

def calculate_center_of_mass(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    atoms = [atom for atom in structure.get_atoms() if atom.element != 'H']
    coordinates = np.array([atom.coord for atom in atoms])
    center_of_mass = np.mean(coordinates, axis=0)
    return center_of_mass

def run_vina(receptor_pdbqt, ligand_pdbqt, out_pdbqt, center_x, center_y, center_z, size_x, size_y, size_z):
    vina_command = f'vina --receptor {receptor_pdbqt} --ligand {ligand_pdbqt} --out {out_pdbqt} --center_x {center_x} --center_y {center_y} --center_z {center_z} --size_x {size_x} --size_y {size_y} --size_z {size_z}'
    os.system(vina_command)
    print(f'Vina 도킹이 완료되었습니다. 결과는 {out_pdbqt} 파일에 저장되었습니다.')

def extract_ligand_center_and_energy(pdbqt_file):
    with open(pdbqt_file, 'r') as file:
        lines = file.readlines()
    
    x_coords, y_coords, z_coords = [], [], []
    energies = []

    for line in lines:
        if line.startswith("ATOM") and line[13:15].strip() == "C":
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            x_coords.append(x)
            y_coords.append(y)
            z_coords.append(z)
        elif line.startswith("REMARK VINA RESULT"):
            parts = line.split()
            energies.append(float(parts[3]))
    
    if not energies:
        raise ValueError("No energies found in the PDBQT file.")
    
    center_x = np.mean(x_coords)
    center_y = np.mean(y_coords)
    center_z = np.mean(z_coords)
    energy = np.mean(energies) if energies else None
    
    return [center_x, center_y, center_z], energy

def main(protein_pdbqt, ligand_pdbqt, protein_pdb):
    # Calculate the center of mass for the protein
    center_x, center_y, center_z = calculate_center_of_mass(protein_pdb)
    print(f'그리드 박스 중심: X={center_x}, Y={center_y}, Z={center_z}')

    # Grid box size (example values, adjust as needed)
    size_x, size_y, size_z = 20, 20, 20

    # Define the grid dimensions (e.g., 4x4x4)
    num_points_per_axis = 5
    grid_points_x = np.linspace(center_x - 40, center_x + 40, num_points_per_axis)
    grid_points_y = np.linspace(center_y - 40, center_y + 40, num_points_per_axis)
    grid_points_z = np.linspace(center_z - 40, center_z + 40, num_points_per_axis)

    ligand_centers = []
    binding_energies = []

    # Random docking runs within each grid volume
    run_index = 0
    for gx in grid_points_x:
        for gy in grid_points_y:
            for gz in grid_points_z:
                adjusted_center_x = gx + np.random.uniform(-10, 10)
                adjusted_center_y = gy + np.random.uniform(-10, 10)
                adjusted_center_z = gz + np.random.uniform(-10, 10)
                
                out_pdbqt = f"docking_output_{run_index+1}.pdbqt"
                
                print(f'Running Vina with center: X={adjusted_center_x}, Y={adjusted_center_y}, Z={adjusted_center_z}')
                run_vina(protein_pdbqt, ligand_pdbqt, out_pdbqt, adjusted_center_x, adjusted_center_y, adjusted_center_z, size_x, size_y, size_z)

                ligand_center, energy = extract_ligand_center_and_energy(out_pdbqt)
                if energy is not None:
                    ligand_centers.append(ligand_center)
                    binding_energies.append(energy)

                run_index += 1

    # Convert to numpy array for clustering and analysis
    ligand_centers = np.array(ligand_centers)
    binding_energies = np.array(binding_energies)

    # Save ligand centers and binding energies for further analysis
    np.savetxt("ligand_centers.csv", ligand_centers, delimiter=",")
    np.savetxt("binding_energies.csv", binding_energies, delimiter=",")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("사용 방법: python random_docking.py <protein_pdbqt> <ligand_pdbqt> <protein_pdb>")
        sys.exit(1)

    protein_pdbqt = sys.argv[1]
    ligand_pdbqt = sys.argv[2]
    protein_pdb = sys.argv[3]
    main(protein_pdbqt, ligand_pdbqt, protein_pdb)
