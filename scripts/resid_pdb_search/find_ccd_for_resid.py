import os
import pandas as pd
import sys
import inspect

# Import config.py
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

import config

def find_common_ccd(csv_file):
    output = pd.read_csv(csv_file)

    # Split values in each cell based on '/'
    output['split_values'] = output['CCD_ID'].str.split('/')

    # Count occurrences of each value
    value_counts = {}
    for row_values in output['split_values']:
        for value in row_values:
            if value not in ["ASN", "MSE", "SER", "THR"]:  # Exclude "ASN" and "MSE"
                value_counts[value] = value_counts.get(value, 0) + 1

    # Calculate percentage of occurrence for each value
    total_rows = len(output)
    threshold = 0.75 * total_rows
    common_values = [value for value, count in value_counts.items() if count > threshold]

    return common_values

def common_values_to_csv(output_file_path):
    directory = os.path.join(config.base_directory, 'resid_pdb_search/outputs/resid_ccd_outputs')
    csv_files = [f for f in os.listdir(directory) if f.endswith(".csv")]

    common_values_list = []
    for csv_file in csv_files:
        file_path = os.path.join(directory, csv_file)
        output = find_common_ccd(file_path)
        filename = os.path.splitext(csv_file)[0]  # Extract filename without extension
        common_values_list.append({'RESID': filename, 'CCD_ID': ','.join(output)})

    # Create DataFrame from list of dictionaries
    df_common_values = pd.DataFrame(common_values_list)
    df_common_values.to_csv(output_file_path, index=False)

    print("Common values saved to common_values.csv")

if __name__ == "__main__":
    
    output_file_path = os.path.join(config.base_directory, 'resid_pdb_search/outputs/extracted_ccd_from_resid.csv')
    common_values_to_csv(output_file_path)
    