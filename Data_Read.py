import os
import re

def extract_info_from_log(file_path):
    scf_done_matches = []
    input_orientation_matches = []

    with open(file_path, 'r') as file:
        content = file.read()

        scf_done_pattern = re.compile(r'SCF Done: ([^\n]+)', re.IGNORECASE)
        input_orientation_pattern = re.compile(r'Input orientation:(.*?)-{70}', re.DOTALL)

        scf_done_matches = scf_done_pattern.findall(content)
        input_orientation_matches_raw = input_orientation_pattern.findall(content)

        for match in input_orientation_matches_raw:
            # Use findall to capture both rows of column names
            column_names_match = re.findall(r'\s+Center\s+Atomic\s+Atomic\s+Coordinates \(Angstroms\)\n\s+Number\s+Number\s+Type\s+X\s+Y\s+Z', match)
            # Use findall to capture the data rows
            input_orientation_match = re.findall(r'\s+(\d+)\s+(\d+)\s+(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)', match)

            if column_names_match and input_orientation_match:
                # Convert coordinates to tuples with the first three columns
                columns_info = [(int(row[0]), int(row[1]), int(row[2])) for row in input_orientation_match]
                coordinates = [(float(row[3]), float(row[4]), float(row[5])) for row in input_orientation_match]
                input_orientation_matches.append((column_names_match, columns_info, coordinates))

    return scf_done_matches, input_orientation_matches

def process_folder(folder_name, output_file_path):
    with open(output_file_path, 'w') as output_file:
        for root, dirs, files in os.walk(folder_name):
            for file in files:
                if file.endswith(".log"):
                    file_path = os.path.join(root, file)
                    file_name = os.path.splitext(file)[0]  # Get filename without extension
                    scf_done_info, input_orientation_info_list = extract_info_from_log(file_path)

                    output_file.write(f"File: {file_name}\n")

                    for i, scf_done_value in enumerate(scf_done_info):
                        output_file.write(f"Iteration {i + 1} - SCF Done Information:\n")
                        output_file.write(f"{scf_done_value}\n")

                        output_file.write("Input Orientation Information:")
                        input_orientation_info = input_orientation_info_list[i] if i < len(input_orientation_info_list) else []

                        if input_orientation_info:
                            # Write the column names (with two rows for formatting)
                            output_file.write("".join(input_orientation_info[0]) + "\n")

                            # Write the first three columns
                            for columns_info, coord in zip(input_orientation_info[1], input_orientation_info[2]):
                                output_file.write("{:<12} {:<12} {:<12} {:<12} {:<12} {:<12}\n".format(*columns_info, *coord))
                            output_file.write('-----------------')
                            output_file.write("\n")
                            
folder_name = 'bdpolarity-KD/BDP-ph'
conclusions_log_path = 'ExtractedData.log'
process_folder(folder_name, conclusions_log_path)