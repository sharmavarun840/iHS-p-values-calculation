import pandas as pd
import numpy as np
from scipy.stats import norm
import glob
import os

# Find any .ihs.out file in the current directory
file_list = glob.glob("*.out")

if not file_list:
    raise FileNotFoundError("No .out file found in the current directory.")

# Use the first .ihs.out file found
file_path = file_list[0]
print(f"Processing file: {file_path}")

# Load the .ihs.out file
data = pd.read_csv(file_path, delim_whitespace=True, header=None)

# Assign column names based on the file structure
data.columns = ['chromosome', 'snp_id', 'position', 'value1', 'value2', 'value3', 'ihs']

# Convert the 'ihs' column to numeric values, forcing any errors to NaN
data['ihs'] = pd.to_numeric(data['ihs'], errors='coerce')

# Drop rows with NaN values in 'ihs' column
data = data.dropna(subset=['ihs'])

# Calculate two-sided p-values from the iHS scores
data['p_value'] = 2 * norm.cdf(-np.abs(data['ihs']))

# Rank the p-values
data = data.sort_values('p_value').reset_index(drop=True)
data['rank'] = data.index + 1

# Adjust p-values for multiple testing using the Benjamini-Hochberg (BH) method
total_tests = len(data)
data['bh_adjusted_p'] = (data['p_value'] * total_tests / data['rank']).clip(upper=1.0)

# Generate the output file name
output_path = os.path.splitext(file_path)[0] + "_with_p_values.out"

# Save the results to a new file
data.to_csv(output_path, sep='\t', index=False)

print(f"P-values and adjusted p-values have been calculated and saved to {output_path}")
