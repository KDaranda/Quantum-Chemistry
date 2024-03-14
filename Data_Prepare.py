import pandas as pd
from decimal import Decimal

def parse_internal_coords(file_path):
    data = []
    molecule_data = {'molecule_name': None, 'definitions': [], 'values': [], 'SCF_energy': None, 'iteration': None}

    with open(file_path, 'r') as file:
        lines = file.readlines()
        for i, line in enumerate(lines):
            if line.startswith('File:'):
                # Capture the molecule name for subsequent iterations
                current_molecule_name = line.split(':')[1].strip()

            elif line.startswith('Iteration:'):
                # Start a new molecule_data dictionary for each iteration
                if molecule_data['definitions']:  # Check if there's already data collected to avoid appending empty data at the start
                    data.append(molecule_data)
                    molecule_data = {'molecule_name': current_molecule_name, 'definitions': [], 'values': [], 'SCF_energy': None, 'iteration': None}  # Reset for new molecule data
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

file_path = 'ConvertedCoordinates.log'
parsed_data = parse_internal_coords(file_path)

print('\n--------------\nDATA IS PARSED\n--------------\n')

# DATA FLATTENING TO A LIST OF TUPLES
def flatten_data(parsed_data):
    flattened_data = []
    for entry in parsed_data:
        molecule_name = entry['molecule_name']
        scf_energy = entry['SCF_energy']
        iteration = entry['iteration']
        for definition, value in zip(entry['definitions'], entry['values']):
            flattened_data.append((molecule_name, definition, value, scf_energy, iteration))
    return flattened_data

# CREATING FEATURE DATAFRAME
def encode_features(flattened_data):
    # Convert flattened data to DataFrame
    columns = ['molecule_name', 'definition', 'value', 'SCF_energy', 'iteration']
    df = pd.DataFrame(flattened_data, columns=columns)

    # Pivot the table to get unique definitions as columns
    feature_df = df.pivot_table(index=['molecule_name', 'SCF_energy', 'iteration'], 
                                columns='definition', 
                                values='value', 
                                fill_value=0).reset_index()

    # Reset index might be necessary if you want 'molecule_name' as a column
    feature_df.reset_index(drop=True, inplace=True)

    return feature_df

# Assuming `parsed_data` is your initial list of dictionaries after parsing the log file
flattened_data = flatten_data(parsed_data)

print('\n-----------------\nDATA IS FLATTENED\n-----------------\n')

feature_df = encode_features(flattened_data)

print(feature_df.head())
'''print('\n--------------------\nDATAFRAME IS CREATED\n--------------------\n')

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split, GridSearchCV
import xgboost as xgb
from sklearn.metrics import mean_squared_error
import numpy as np

# Define features and target
X = feature_df.drop(['SCF_energy', 'molecule_name', 'iteration'], axis=1)
y = feature_df['SCF_energy']

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

print('\n-------------\nDATA IS SPLIT\n-------------\n')

# Normalize feature data
scaler_X = StandardScaler()

X_train_scaled = scaler_X.fit_transform(X_train)
X_test_scaled = scaler_X.transform(X_test)

print('\n------------------\nDATA IS NORMALIZED\n------------------\n')

# Initialize XGBRegressor
xgb_reg = xgb.XGBRegressor(objective='reg:squarederror')

# Setup the hyperparameter grid for GridSearchCV
param_grid = {
    'max_depth': [3, 30, 300],
    'eta': [0.05, 0.1, 0.3],
    'subsample': [0.7, 0.8, 0.9],
    'colsample_bytree': [0.7, 0.8, 0.9],
}

# Initialize GridSearchCV
grid_search = GridSearchCV(estimator=xgb_reg, param_grid=param_grid, cv=5, scoring='neg_root_mean_squared_error', verbose=3)

print('\n---------------------\nSTARTING THE TRAINING\n---------------------\n')

# Fit GridSearchCV
grid_search.fit(X_train_scaled, y_train)

# Best parameters
print("Best parameters found: ", grid_search.best_params_)

# Best model
best_model = grid_search.best_estimator_

# Predictions
predictions = best_model.predict(X_test_scaled)

# Calculate RMSE
rmse = np.sqrt(mean_squared_error(y_test, predictions))
print(f"RMSE: {rmse}")'''

'''
# Convert data to DMatrix format for XGBoost
dtrain = xgb.DMatrix(X_train_scaled, label=y_train) 
dtest = xgb.DMatrix(X_test_scaled, label=y_test) 

# Parameters for XGBoost
params = {
    'max_depth': 100,
    'eta': 0.1,
    'objective': 'reg:squarederror',
    'eval_metric': 'rmse'
}

print('\n---------------------\nSTARTING THE TRAINING\n---------------------\n')

num_rounds = 1000
evals_result = {}
bst = xgb.train(params, dtrain, num_rounds, evals=[(dtrain, 'train'), (dval, 'validation')], early_stopping_rounds=10, evals_result = evals_result, verbose_eval=True)

# Predictions
predictions = bst.predict(dtest)

# Calculate and print RMSE
rmse = np.sqrt(mean_squared_error(y_test.astype(float), predictions))  # Ensure y_test is float
print(f"RMSE: {rmse}")

# Optionally, print evaluation metrics
print('\nFinal testing metrics:')
print(evals_result['validation']['rmse'][-1])  # Example: printing final RMSE for test set
'''