# visualize_energy_pymol.py

import sys
import numpy as np

def extract_ligand_centers_and_energies(pdbqt_files):
    ligand_centers = []
    binding_energies = []
    for pdbqt_file in pdbqt_files:
        with open(pdbqt_file, 'r') as file:
            lines = file.readlines()
        
        x_coords, y_coords, z_coords = [], [], []
        energy = None
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
                energy = float(parts[3])
        
        if x_coords and y_coords and z_coords:
            center_x = np.mean(x_coords)
            center_y = np.mean(y_coords)
            center_z = np.mean(z_coords)
            ligand_centers.append([center_x, center_y, center_z])
            binding_energies.append(energy)
    
    return np.array(ligand_centers), np.array(binding_energies)

def assign_colors(binding_energies):
    sorted_indices = np.argsort(binding_energies)
    num_ligands = len(binding_energies)
    
    colors = np.zeros((num_ligands, 3), dtype=int)
    color_map = [
        [255, 0, 0],    # Red
        [255, 165, 0],  # Orange
        [255, 255, 0],  # Yellow
        [0, 128, 0],    # Green
        [0, 0, 255]     # Blue
    ]

    # Divide the ligands into 5 groups based on binding energy rank
    for rank, index in enumerate(sorted_indices):
        color_index = rank * len(color_map) // num_ligands
        colors[index] = color_map[color_index]

    return colors

def save_pymol_script(protein_file, ligand_centers, colors):
    with open("ligand_energies.pml", "w") as f:
        f.write(f"load {protein_file}\n")
        for i, (center, color) in enumerate(zip(ligand_centers, colors)):
            f.write(f"pseudoatom ligand_center_{i+1}, pos=[{center[0]},{center[1]},{center[2]}], color=[{color[0]},{color[1]},{color[2]}]\n")
        f.write("show spheres, ligand_center_*\n")
        f.write("set sphere_scale, 1.0, ligand_center_*\n")  # Increase sphere scale for larger visualization
        for i in range(len(ligand_centers)):
            f.write(f"set_color color_{i+1}, [{colors[i][0]/255.0},{colors[i][1]/255.0},{colors[i][2]/255.0}]\n")
            f.write(f"color color_{i+1}, ligand_center_{i+1}\n")

def main(protein_file, pdbqt_files):
    ligand_centers, binding_energies = extract_ligand_centers_and_energies(pdbqt_files)
    
    colors = assign_colors(binding_energies)
    save_pymol_script(protein_file, ligand_centers, colors)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("사용 방법: python visualize_energy_pymol.py <단백질 PDB 파일 경로> <PDBQT 파일 경로>...")
        sys.exit(1)

    protein_file = sys.argv[1]
    pdbqt_files = sys.argv[2:]
    main(protein_file, pdbqt_files)
