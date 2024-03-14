import math_connect

def extract_cartesian_coordinates(log_file_path):
    coordinates_list = []
    file_info = ""  # To capture the "File:" line
    scf_energy = None 
    coordinates = []
    start_reading = False

    with open(log_file_path, 'r') as file:
        for line in file:
            if line.startswith('File:'):
                file_info = line.strip()  # Capture the "File:" line
                #scf_energy = None

            if line.startswith(" E(RM062X) =  "):  # Example pattern, adjust if needed
                scf_energy = line.strip()[:32]  # Extracting SCF energy value

            if 'Input Orientation Information:' in line:
                start_reading = True
                # If starting a new section, reset the current coordinates and save the file info
                if coordinates:
                    coordinates_list.append((file_info, scf_energy, coordinates))
                    coordinates = []
                continue

            if start_reading and '-----------------' in line:
                start_reading = False
                continue

            if start_reading:
                parts = line.split()
                try:
                    x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
                    coordinates.append((x, y, z))
                except (IndexError, ValueError):
                    continue

    if coordinates:
        coordinates_list.append((file_info, scf_energy, coordinates))

    return coordinates_list

def convert_cartesian_to_internal(all_cartesian_coords):
    all_internal_coords = []
    
    for file_info, scf_energy, cartesian_coords in all_cartesian_coords:
        coords = {index: (x, y, z) for index, (x, y, z) in enumerate(cartesian_coords, start=1)}

        bonds = math_connect.get_bonds(coords)
        angles = math_connect.get_angles(bonds, coords)
        dihedrals = math_connect.get_diheds(angles, coords)

        internal_coords = []

        for bond, length in bonds.items():
            internal_coords.append(('R', bond, length))

        for angle, value in angles.items():
            internal_coords.append(('A', angle, value))

        for dihedral, value in dihedrals.items():
            internal_coords.append(('D', dihedral, value))
        
        all_internal_coords.append((file_info, scf_energy, internal_coords))
    
    return all_internal_coords

def write_internal_coordinates(all_internal_coords, output_file):
    previous_file = " "
    it_count = 0
    with open(output_file, 'w') as file:
        for file_info, scf_energy, internal_coords in all_internal_coords:
            if previous_file != file_info:
                previous_file = file_info
                it_count = 1
                file.write(f"{file_info}\n")
            else: it_count += 1
            file.write(f"Iteration: {it_count}\n")
            file.write(f"{scf_energy}\n")
            file.write("-----------------------------------------------\n")
            file.write("! Name  Definition              Value         !\n")
            file.write("-----------------------------------------------\n")
            
            r_counter, a_counter, d_counter = 1, 1, 1

            # Separating and sorting internal coordinates based on their types and indices
            r_coords = sorted([x for x in internal_coords if x[0] == 'R'], key=lambda x: x[1])
            a_coords = sorted([x for x in internal_coords if x[0] == 'A'], key=lambda x: x[1])
            d_coords = sorted([x for x in internal_coords if x[0] == 'D'], key=lambda x: x[1])

            for _, indices, value in r_coords:
                definition = f"R({','.join(map(str, indices))})"
                file.write(f"! R{r_counter:<5} {definition:<20} {value:<15.6f} !\n")
                r_counter += 1

            for _, indices, value in a_coords:
                definition = f"A({','.join(map(str, indices))})"
                file.write(f"! A{a_counter:<5} {definition:<20} {value:<15.6f} !\n")
                a_counter += 1

            for _, indices, value in d_coords:
                definition = f"D({','.join(map(str, indices))})"
                file.write(f"! D{d_counter:<5} {definition:<20} {value:<15.6f} !\n")
                d_counter += 1

            file.write("-----------------------------------------------\n")
            file.write("\n")


# Main execution
log_file_path = 'ExtractedData.log'
converted_file_path = 'ConvertedCoordinates.log'

# Extracting Cartesian coordinates
cartesian_coordinates = extract_cartesian_coordinates(log_file_path)

# Converting to internal coordinates
internal_coordinates = convert_cartesian_to_internal(cartesian_coordinates)

# Writing converted coordinates to a log file
write_internal_coordinates(internal_coordinates, converted_file_path)