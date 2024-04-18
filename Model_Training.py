from Data_Prepare import filtered_df, specific_molecule
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split, GridSearchCV
import xgboost as xgb
from sklearn.metrics import mean_absolute_error
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

scaler_X = StandardScaler()

# Parameters for XGBoost booster model
params = {
    'max_depth': 5,
    'eta': 0.3,
    'subsample': 0.5,
    'objective': 'reg:squarederror',
    'eval_metric': 'mae',
    'colsample_bytree': 0.1,
    'alpha': 0,
    'gamma': 0,
    'lambda': 2,
    'min_child_weight': 2,
    #'n_estimators': 300
}

molecule_solvent_pairs = filtered_df[['molecule', 'solvent']].drop_duplicates()
pair_count = len(molecule_solvent_pairs)

# Determine the layout of the subplots
cols = 3  # Number of columns in the subplot grid
rows = (pair_count + cols - 1) // cols  # Calculate rows needed

fig, axs = plt.subplots(rows, cols, figsize=(15, rows * 5))
axs = axs.flatten()  # Flatten to 1D array for easy iteration

mae_total = 0
mae_count = 0

for index, (ax, (_, row)) in enumerate(zip(axs, molecule_solvent_pairs.iterrows())):
    molecule = row['molecule']
    solvent = row['solvent']

    # Filter data for the current molecule-solvent pair
    df_pair = filtered_df[(filtered_df['molecule'] == molecule) & (filtered_df['solvent'] == solvent)]
    X = df_pair.drop(['SCF_energy', 'molecule', 'solvent', 'iteration', 'energy_diff'], axis=1)
    y = df_pair['energy_diff']

    # Split the data into training and testing sets
    X_train, X_ham, y_train, y_ham = train_test_split(X, y, test_size=0.4, random_state=42)
    X_test, X_val, y_test, y_val = train_test_split(X_ham, y_ham, test_size=0.5, random_state=42)

    print('\n-------------\nDATA IS SPLIT\n-------------\n')

    X_train_scaled = scaler_X.fit_transform(X_train)
    X_test_scaled = scaler_X.transform(X_test)
    X_val_scaled = scaler_X.transform(X_val)

    print('\n------------------\nDATA IS NORMALIZED\n------------------\n')

    # Convert data to DMatrix format for XGBoost
    dtrain = xgb.DMatrix(X_train_scaled, label=y_train) 
    dtest = xgb.DMatrix(X_test_scaled, label=y_test)
    dval = xgb.DMatrix(X_val_scaled, label=y_val)

    print('\n---------------------\nSTARTING THE TRAINING\n---------------------\n')

    # Train the model
    evals_result = {}
    bst = xgb.train(params, dtrain, 100, evals=[(dtrain, 'train'), (dtest, 'test'), (dval, 'validation')], early_stopping_rounds=5, evals_result=evals_result, verbose_eval=True)

    # Predictions
    predictions = bst.predict(dtest)

    # Calculate MAE
    mae = mean_absolute_error(y_test.astype(float), predictions)

    # Optionally, print evaluation metrics
    print('\nFinal best validation subset MAE for ' + molecule + "-" + solvent)
    print(min(evals_result['validation']['mae']))  # Example: printing final MAE for test set

    print('\nFinal best actual MAE for ' + molecule + "-" + solvent)
    print(min(evals_result['test']['mae']))  # Example: printing final MAE for test set
    
    # Store mae for future calculations
    mae_total = mae_total + mae
    mae_count = mae_count + 1

    # Plotting the training and test metrics
    epochs = len(evals_result['train']['mae'])
    x_axis = range(0, epochs)
    train_mae = evals_result['train']['mae']
    test_mae = evals_result['test']['mae']
    validation_mae = evals_result['validation']['mae']

    ax.plot(x_axis, evals_result['train']['mae'], label='Train', color = 'red')
    ax.plot(x_axis, evals_result['test']['mae'], label='Test', color = 'blue')
    ax.plot(x_axis, evals_result['validation']['mae'], label='Validation', color = 'green')

    # Highlight the minimum points on the graph
    min_train_idx = np.argmin(train_mae)
    min_test_idx = np.argmin(test_mae)
    min_validation_idx = np.argmin(validation_mae)

    ax.scatter(min_train_idx, train_mae[min_train_idx], color='blue', s=20, label=f'Train Min = {train_mae[min_train_idx]:.5f}')
    ax.scatter(min_test_idx, test_mae[min_test_idx], color='red', s=20, label=f'Test Min = {test_mae[min_test_idx]:.5f}')
    ax.scatter(min_validation_idx, validation_mae[min_validation_idx], color='black', s=20, label=f'Validation Min = {validation_mae[min_validation_idx]:.5f}')

    ax.set_title(f"{molecule}-{solvent}")
    ax.set_xlabel('Epochs')
    ax.set_ylabel('MAE')
    ax.legend()


# Adjust layout
plt.tight_layout()
plt.show()

average_mae = mae_total/mae_count
print(f"\nAverage MAE for {specific_molecule} with different solvents is equal to {average_mae}")


# Define features and target
'''X = filtered_df.drop(['SCF_energy', 'molecule', 'solvent', 'iteration', 'energy_diff'], axis=1)
y = filtered_df['energy_diff']

# Split the data into training and testing sets
X_train, X_ham, y_train, y_ham = train_test_split(X, y, test_size=0.4, random_state=42)
X_test, X_val, y_test, y_val = train_test_split(X_ham, y_ham, test_size=0.5, random_state=42)

print('\n-------------\nDATA IS SPLIT\n-------------\n')

# Normalize feature data
X_train_scaled = scaler_X.fit_transform(X_train)
X_test_scaled = scaler_X.transform(X_test)
X_val_scaled = scaler_X.transform(X_val)

print('\n------------------\nDATA IS NORMALIZED\n------------------\n')

# Convert data to DMatrix format for XGBoost
dtrain = xgb.DMatrix(X_train_scaled, label=y_train) 
dtest = xgb.DMatrix(X_test_scaled, label=y_test) 
dval = xgb.DMatrix(X_val_scaled, label=y_val)

print('\n---------------------\nSTARTING THE TRAINING\n---------------------\n')

evals_result = {}
bst = xgb.train(params, dtrain, num_rounds = 100, evals=[(dtrain, 'train'), (dval, 'validation')], early_stopping_rounds=5, evals_result = evals_result, verbose_eval=True)

# Predictions
predictions = bst.predict(dtest)

# Calculate and print MAE
mae = mean_absolute_error(y_test.astype(float), predictions)  # Ensure y_test is float
print(f"MAE for testing subset: {mae}")

# Optionally, print evaluation metrics
print('\nFinal testing metrics for ' + specific_molecule + ' in validation subset')
print(evals_result['validation']['mae'][-6])  # Example: printing final MAE for test set

print('\n---------\nPLOTTING\n---------\n')

# Plotting the training and validation metrics
epochs = len(evals_result['train']['mae'])
x_axis = range(0, epochs)
fig, ax = plt.subplots()
ax.plot(x_axis, evals_result['train']['mae'], label='Train')
ax.plot(x_axis, evals_result['validation']['mae'], label='Validation')

# Highlight the minimum MAE for both training and validation
min_train_idx = np.argmin(evals_result['train']['mae'])
min_val_idx = np.argmin(evals_result['validation']['mae'])
ax.scatter(min_train_idx, evals_result['train']['mae'][min_train_idx], color='red', label='Train Min', zorder=5)
ax.scatter(min_val_idx, evals_result['validation']['mae'][min_val_idx], color='blue', label='Validation Min', zorder=5)

# Annotate the minimum points
ax.annotate(f'{evals_result["train"]["mae"][min_train_idx]:.4f}', (min_train_idx, evals_result['train']['mae'][min_train_idx]), textcoords="offset points", xytext=(-15,-10), ha='center', color='blue')
ax.annotate(f'{evals_result["validation"]["mae"][min_val_idx]:.4f}', (min_val_idx, evals_result['validation']['mae'][min_val_idx]), textcoords="offset points", xytext=(-15,10), ha='center', color='red')

ax.legend()
plt.ylabel('MAE')
plt.title('XGBoost MAE')
plt.text(0, 0.01, f"MAE for testing subset: {mae}")
plt.show()'''


'''# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.4, random_state=42)

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
    'max_depth': [5],
    'min_child_weight': [2],
    'subsample': [0.3],
    'colsample_bytree': [0.1],
    'eta': [0.1],
    'lambda': [2],
    'alpha': [0],
    'gamma': [0],
    'n_estimators': [300]
}


# Initialize GridSearchCV
grid_search = GridSearchCV(estimator=xgb_reg, param_grid=param_grid, cv=5, scoring='neg_mean_absolute_error', verbose=3, n_jobs=-1, error_score='raise')
print('\n---------------------\nSTARTING THE TRAINING\n---------------------\n')

# Fit GridSearchCV
grid_search.fit(X_train_scaled, y_train)

# Best parameters
print("Best parameters found: ", grid_search.best_params_)

# Best model
best_model = grid_search.best_estimator_
#best_model.save_model(0001.model)

# Predictions
predictions = best_model.predict(X_test_scaled)

# Calculate RMSE
mae = mean_absolute_error(y_test, predictions)
print(f"MAE for testing subset: {mae}")'''
