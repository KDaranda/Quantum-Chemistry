import h5py
import numpy as np

def parse_cartesian_file(file_path):
    energies = []
    coordinates = []
    atom_types = []
    names = []
    current_coords = []
    current_atoms = []
    current_name = []

    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if line.startswith('File:'):
                current_name = '-'.join(line.split()[1].split('-')[:3])  # Remove part after the third dash
            elif line.startswith(' E(RM062X) ='):
                energy = float(line.split('=')[1].strip().split()[0])
                energies.append(energy)
                names.append(current_name)
                if current_coords:
                    coordinates.append(np.array(current_coords))
                    atom_types.append(np.array(current_atoms))
                current_coords = []
                current_atoms = []
            elif line.lstrip().startswith(' Center '):
                continue
            elif line.lstrip().startswith(' Number '):
                continue
            elif any(char.isdigit() for char in line.split()[0]):  # Simple check for coordinate line
                parts = line.split()
                atom_type = int(parts[1])  # Assuming second part is atom number (atomic number)
                x, y, z = map(float, parts[-3:])
                current_atoms.append(atom_type)
                current_coords.append([x, y, z])
        
        if current_coords:
            coordinates.append(np.array(current_coords))
            atom_types.append(np.array(current_atoms))

    return names, atom_types, coordinates, energies

def create_hdf5_dataset(names, atom_types, coordinates, energies, output_file):
    with h5py.File(output_file, 'w') as f:
        grp = f.create_group('molecules')
        for i, (name, coords, atoms, energy) in enumerate(zip(names, coordinates, atom_types, energies)):
            mol_grp = grp.create_group(f'molecule_{i}')
            mol_grp.create_dataset('name', data=name)
            mol_grp.create_dataset('atoms', data=atoms)
            mol_grp.create_dataset('coordinates', data=coords)
            mol_grp.create_dataset('energy', data=energy)

# Parse the Cartesian coordinates file
names, atom_types, coordinates, energies = parse_cartesian_file('ExtractedData.log')

# Create the HDF5 dataset for ANI-2x
hdf5_output_path = 'ani2x_data_with_names.h5'
create_hdf5_dataset(names, atom_types, coordinates, energies, hdf5_output_path)

# Confirmation message
print(hdf5_output_path, len(energies))

'''SKAITOM DATASET'''

def read_hdf5_file(file_path, num_entries):
    with h5py.File(file_path, 'r') as f:
        entries = []
        for i, key in enumerate(f['molecules'].keys()):
            if i >= num_entries:
                break
            mol_grp = f['molecules'][key]
            name = mol_grp['name'][()].decode('utf-8')
            atoms = mol_grp['atoms'][:]
            coordinates = mol_grp['coordinates'][:]
            energy = mol_grp['energy'][()]
            entries.append((name, atoms, coordinates, energy))
        return entries

# Reading the first 100 entries from the HDF5 file
first_100_entries = read_hdf5_file(hdf5_output_path, 5)

# Displaying the first 3 entries for brevity
print(first_100_entries[:3])
