import os
import csv
import pandas as pd
import re
import sys
import inspect

# Import config.py
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

import config

# Parse uniprot ptmlist.txt
input_file_path = os.path.join(config.base_directory, 'uniprot_ptm_list/inputs/ptmlist.txt')
output_csv_path = os.path.join(config.base_directory, 'uniprot_ptm_list/outputs/raw_uniprot_ptm_list.csv')

with open(input_file_path, 'r') as input_file, open(output_csv_path, 'w', newline='') as output_csv:

    column_order = ['ID', 'AC', 'FT', 'TG', 'PA', 'PP', 'CF', 'MM', 'MA', 'LC', 'TR', 'KW', 'DR', 'ChEBI', 'RESID']
    csv_writer = csv.DictWriter(output_csv, fieldnames=column_order)
    csv_writer.writeheader()

    entry_data = {}

    for _ in range(47): # Skip the first 47 blank lines
        next(input_file)

    # Define a function to extract ChEBI number
    def extract_chebi_number(dr_string):
        match = re.search(r'CHEBI:(\d+)', dr_string)
        if match:
            return match.group(1)
        else:
            return None

    def extract_resid_number(dr_string):
        match = re.search(r'RESID; (\w+)', dr_string)
        if match:
            return match.group(1)
        else:
            return None

    for line in input_file:
        line = line.strip()

        if not line:
            continue

        # Check if the line ends with "//" which is the end of a block
        if line.endswith('//'):
            entry_data['ChEBI'] = ', '.join(filter(None, map(extract_chebi_number, entry_data.get('DR', []))))
            entry_data['RESID'] = ', '.join(filter(None, map(extract_resid_number, entry_data.get('DR', []))))

            csv_writer.writerow(entry_data)
            entry_data = {}
        else:
            if ' ' in line:
                key, value = line.split(maxsplit=1)
                entry_data.setdefault(key, []).append(value)

print(f'CSV file successfully created at: {output_csv_path}')

input_file_path = os.path.join(config.base_directory, 'uniprot_ptm_list/outputs/raw_uniprot_ptm_list.csv')
output_file_path = os.path.join(config.base_directory, 'uniprot_ptm_list/outputs/final_uniprot_ptm_list.csv')

output_uniprot = pd.read_csv(input_file_path)

output_uniprot = output_uniprot.astype(str)
output_uniprot.replace('', pd.NA, inplace=True)

columns_to_strip = output_uniprot.columns[:12]  # Adjust this range if needed
for column in columns_to_strip:
    output_uniprot[column] = output_uniprot[column].str.replace("[\[\]]", "", regex=True)

for col in output_uniprot.columns[:12]:
    output_uniprot[col] = output_uniprot[col].str.replace("'", "")

column_names = ["ID", "target_residue", "residue_location", "protein_location", "PTM_type", "ChEBI", "RESID"]
final_ptm_list = output_uniprot[["ID", "TG", "PA", "PP", "KW", "ChEBI", "RESID"]].copy()
final_ptm_list.columns = column_names

final_ptm_list.to_csv(output_file_path, index=False)