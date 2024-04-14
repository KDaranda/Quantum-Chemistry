from decimal import Decimal
import pandas as pd
import polars as pl 

def parse_internal_coords(file_path):
    data = []
    molecule_data = {'molecule': None, 'solvent': None, 'definitions': [], 'values': [], 'SCF_energy': None, 'iteration': None}

    with open(file_path, 'r') as file:
        lines = file.readlines()
        for i, line in enumerate(lines):
            if line.startswith('File:'):
                file_name = line.split(':')[1].strip() 
                parts = file_name.split('-') # Assuming naming convention is 'molecule-solvent-...'
                molecule, solvent = parts[0] + "-" + parts[1], parts[2]
            
            elif line.startswith('Iteration:'):
                # Start a new molecule_data dictionary for each iteration
                if molecule_data['definitions']:  # Check if there's already data collected to avoid appending empty data at the start
                    data.append(molecule_data)
                    molecule_data = {'molecule': molecule, 'solvent': solvent, 'definitions': [], 'values': [], 'SCF_energy': None, 'iteration': None}  # Reset for new molecule data
                molecule_data['iteration'] = int(line.split(':')[1].strip())

            elif line.startswith('E(RM062X) ='):
                molecule_data['SCF_energy'] = Decimal(line.split("=")[-1].strip())

            elif line.strip().startswith('! R') or line.strip().startswith('! A') or line.strip().startswith('! D'):
                parts = line.strip().split()
                try:
                    # Trying to convert last part to float to ensure it's a value line
                    value = float(parts[3])
                    definition = ' '.join(parts[1:3])  # Adjust based on the actual structure
                    molecule_data['definitions'].append(definition)
                    molecule_data['values'].append(value)
                except ValueError:
                    # Skipping lines that don't contain convertible numeric values
                    continue

    # Adding the last molecule's data
    if molecule_data['definitions']:
        data.append(molecule_data)

    return data

# DATA FLATTENING TO A LIST OF TUPLES
def flatten_data(parsed_data):
    flattened_data = []
    for entry in parsed_data:
        for definition, value in zip(entry['definitions'], entry['values']):
            flattened_data.append((entry['molecule'], entry['solvent'], definition, value, entry['SCF_energy'], entry['iteration']))
    return flattened_data

# CREATING FEATURE DATAFRAME
def create_feature_df(flattened_data):
    # Convert flattened data to DataFrame
    df = pd.DataFrame(flattened_data, columns=['molecule', 'solvent', 'definition', 'value', 'SCF_energy', 'iteration'])

    #Introducing the energy_diff parameter for versatility
    df['min_energy'] = df.groupby(['molecule', 'solvent'])['SCF_energy'].transform('min')
    df['energy_diff'] = df['SCF_energy'] - df['min_energy']

    # Pivot the table to get unique definitions as columns
    feature_df = df.pivot_table(index=['molecule', 'solvent', 'SCF_energy', 'energy_diff', 'iteration'], 
                                columns='definition', 
                                values='value', 
                                fill_value=0).reset_index()

    # Reset index might be necessary if you want 'molecule_name' as a column
    feature_df.reset_index(drop=True, inplace=True)

    filtered_df = feature_df[(feature_df['molecule'] == specific_molecule)]

    return filtered_df



specific_molecule = 'bdp-phenyl'
#pecific_solvent = 'tolu'

file_path = 'ConvertedCoordinates.log'
parsed_data = parse_internal_coords(file_path)

print('\n--------------\nDATA IS PARSED\n--------------\n')

# Assuming `parsed_data` is your initial list of dictionaries after parsing the log file
flattened_data = flatten_data(parsed_data)

print('\n-----------------\nDATA IS FLATTENED\n-----------------\n')

filtered_df = create_feature_df(flattened_data)

print(filtered_df.head())
print('\n--------------------\nDATAFRAME IS CREATED\n--------------------\n')