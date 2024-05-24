import os
import csv
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import numpy as np
import sys
import inspect
from matplotlib.ticker import ScalarFormatter
matplotlib.rcParams['font.family'] = 'Arial'

# Import config.py
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
parentparentdir = os.path.dirname(parentdir)

sys.path.insert(0, parentparentdir) 

import config

def pdb_value_counts_dictionary(directory):
    file_lengths = {}
    for filename in os.listdir(directory):
        if filename != ".DS_Store":  # Skip .DS_Store files
            if os.path.isfile(os.path.join(directory, filename)):
                with open(os.path.join(directory, filename), 'r') as file:
                    reader = csv.reader(file)
                    row_count = sum(1 for row in reader) - 1  # Subtract 1 for the header
                    filename_ccd = os.path.splitext(filename)[0] # remove .csv
                    file_lengths[filename_ccd] = row_count
    
    return file_lengths

def get_smPTM_dict(ptm_ccd_polymer_dict, ptm_ccd_cov_linked_dict, ptm_list, ptm_list_cov_linked):
    
    sm_ptm_list = ptm_list[ptm_list['Lipid'] == 'No']
    sm_ptm_list_cov_linked = ptm_list_cov_linked[ptm_list_cov_linked['Mod_type'] != 'LIPID']

    # Merge ptm_list (in polymer) and ptm_list (cov_linked)
    sm_ptm_list_cov_linked.rename(columns={'CCD': 'CCD_combined'}, inplace=True)
    sm_ptm_list = pd.concat([sm_ptm_list, sm_ptm_list_cov_linked], ignore_index=True)  
    sm_ptm_list['CCD_combined_list'] = sm_ptm_list.groupby('ID')['CCD_combined'].transform(lambda x: ','.join(x)) # If an entry has both in polymer and cov_linked, add it to one value
    sm_ptm_list.drop_duplicates(subset=['ID'], inplace=True)

    PTM_combined_dict = {}
    PTM_combined_dict.update(ptm_ccd_polymer_dict)
    PTM_combined_dict.update(ptm_ccd_cov_linked_dict)

    unique_items = set()

    # Find all CCD values for smPTMs
    for index, row in sm_ptm_list.iterrows():
        # Split the CCD_combined_list string into a list of individual items
        ccd_list = row['CCD_combined_list'].split(', ')
        # Update the set of unique items with the items in ccd_list
        unique_items.update(ccd_list)

    smPTM_dict = {}

    # Create a new dictionary that only contains smPTM CCD values
    for key, value in PTM_combined_dict.items():
        if key in unique_items:
            smPTM_dict[key] = value

    return smPTM_dict

def get_smPTM_type_dict(smPTM_dict):
    
    sm_ptm_list = ptm_list[ptm_list['Lipid'] == 'No']
    sm_ptm_list_cov_linked = ptm_list_cov_linked[ptm_list_cov_linked['Mod_type'] != 'LIPID']

    # Merge ptm_list (in polymer) and ptm_list (cov_linked)
    sm_ptm_list_cov_linked.rename(columns={'CCD': 'CCD_combined'}, inplace=True)
    sm_ptm_list = pd.concat([sm_ptm_list, sm_ptm_list_cov_linked], ignore_index=True)  
    sm_ptm_list['CCD_combined_list'] = sm_ptm_list.groupby('ID')['CCD_combined'].transform(lambda x: ','.join(x)) # If an entry has both in polymer and cov_linked, add it to one value
    sm_ptm_list.drop_duplicates(subset=['ID'], inplace=True)

    sm_ptm_list = sm_ptm_list.astype(str)
    sm_ptm_list['PTM_type'] = sm_ptm_list['PTM_type'].str.rstrip('.')
    sm_ptm_list['PDB_count'] = sm_ptm_list['CCD_combined'].str.split(', ').explode().map(smPTM_dict).groupby(level=0).sum()
    sum_dict = {}

    for keyword in sm_ptm_list['PTM_type'].unique():
        filtered_df = sm_ptm_list[sm_ptm_list['PTM_type'].str.contains(keyword, na=False)]
        pdb_count_sum = filtered_df['PDB_count'].sum()
        sum_dict[keyword] = pdb_count_sum

    sum_dict = {key: value for key, value in sum_dict.items() if key != 'nan'}
        
    # print(sum_dict)

    return sum_dict

def get_lipids_dict(ptm_ccd_polymer_dict, ptm_ccd_cov_linked_dict, ptm_list, ptm_list_cov_linked):

    lipid_polymer_list = ptm_list[ptm_list['Lipid'] == 'Yes']
    lipid_cov_linked_list = ptm_list_cov_linked[ptm_list_cov_linked['Mod_type'] == 'LIPID']
    
    lipid_polymer_ccd_set = set(lipid_polymer_list['CCD_combined'])
    lipid_cov_linked_ccd_set = set(lipid_cov_linked_list['CCD'])

    lipid_polymer_ccd_dict = {key: ptm_ccd_polymer_dict[key] for key in lipid_polymer_ccd_set if key in ptm_ccd_polymer_dict}
    lipid_cov_linked_ccd_dict = {key: ptm_ccd_cov_linked_dict[key] for key in lipid_cov_linked_ccd_set if key in ptm_ccd_cov_linked_dict}
    
    lipid_combined_dict = {}
    lipid_combined_dict.update(lipid_polymer_ccd_dict)
    lipid_combined_dict.update(lipid_cov_linked_ccd_dict)

    lipid_ccd_mapping = {'MYK': "MYR",
                         'P1L': "PLM",
                         'GER': 'FAR'}
    
    merged_dict = {}
    for key, value in lipid_combined_dict.items():
        merged_key = lipid_ccd_mapping.get(key, key)
        merged_dict[merged_key] = merged_dict.get(merged_key, 0) + value

    lipid_keyword_mapping = {'Prenylation': "FAR",
                            'Myristoylation': 'MYR',
                            'Palmitoylation': 'PLM',
                            'Octanoylation': '1IC',
                            'Stearoylation': 'STE',
                            'Cholesterylation': 'CLR',
                            'Palmiteoylation': 'PAM',
                            'S-diacylglycerol': 'DGA'}
    
    lipid_dict = {key: merged_dict[value] for key, value in lipid_keyword_mapping.items()}

    return lipid_dict

def get_glycan_dict():
    glycan_list = pd.read_json(os.path.join(config.base_directory, 'pdb_api_search/inputs/glycosylation_per_year.json'))
    glycan_df = glycan_list.T
    total_nglycans = glycan_df['nGlycans'].sum()
    total_oglycans = glycan_df['oGlycans'].sum()
    total_sglycans = glycan_df['sGlycans'].sum()
    total_cglycans = glycan_df['cGlycans'].sum()

    glycan_stats_dict = {"N-linked": total_nglycans,
                         "O-linked": total_oglycans,
                         "S-linked": total_sglycans,
                         "C-linked": total_cglycans,}
    
    return glycan_stats_dict

def get_total_ptm_dict():
    lipid_dict = get_lipids_dict(ptm_ccd_polymer_dict, ptm_ccd_cov_linked_dict, ptm_list, ptm_list_cov_linked)
    smPTM_type_dict = get_smPTM_type_dict(smPTM_dict)
    glyco_dict = get_glycan_dict()

    combined_dict = {}
    combined_dict.update(lipid_dict)
    combined_dict.update(smPTM_type_dict)
    combined_dict.update(glyco_dict)

    combined_dict['N-glycosylation'] = combined_dict.pop('N-linked')
    combined_dict['O-glycosylation'] = combined_dict.pop('O-linked')
    combined_dict['Phosphorylation'] = combined_dict.pop('Phosphoprotein')

    return combined_dict

def get_removed_entries_dict():
    removed_entries_ccd_polymer_dict = pdb_value_counts_dictionary(removed_entries_ccd_polymer_search_output_directory)
    removed_entires_ccd_cov_linked_dict = pdb_value_counts_dictionary(removed_entries_ccd_cov_linked_search_output_directory)
    removed_entries_ccd_polymer_dict.update(removed_entires_ccd_cov_linked_dict)
    
    return removed_entries_ccd_polymer_dict

def graph(df_dict, title, xaxis_label, yaxis_label):
    
    if "Pyrrolidone carboxylic acid" in df_dict:
        df_dict["Pyroglutamic acid"] = df_dict.pop("Pyrrolidone carboxylic acid")

    top_10_values = dict(sorted(df_dict.items(), key=lambda item: item[1], reverse=True)[:10])
    top_10_values = {key: int(value) for key, value in top_10_values.items()}
    top_10_values = {key: value for key, value in top_10_values.items() if value != 0}

    plt.figure(figsize=(12, 9))  # Increase figure size for better visualization
    bars = plt.bar(top_10_values.keys(), top_10_values.values(), color='#6baed6')

    plt.grid(axis='y', linestyle='-', alpha=0.5, color='gray')
    plt.gca().set_axisbelow(True)  # Move grid lines behind bars
    plt.xlabel(xaxis_label, fontsize=22, labelpad=15)
    plt.ylabel(yaxis_label, fontsize=22)
    plt.xticks(rotation=45, ha='right', fontsize=20)
    plt.yticks(fontsize=20)

    # Show the box around the plot
    plt.gca().spines['top'].set_visible(True)
    plt.gca().spines['right'].set_visible(True)
    plt.gca().spines['left'].set_linewidth(0.5)  # Adjust line width for left axis
    plt.gca().spines['bottom'].set_linewidth(0.5)  # Adjust line width for bottom axis

    plt.tight_layout()  # Adjust layout for better spacing
    plt.savefig(os.path.join(config.base_directory, 'pdb_api_search/analysis/Outputs/Graphs', f'{title}.png'), dpi=600)

    plt.show()

def log_graph(df_dict, title, xaxis_label, yaxis_label):
    top_10_values = dict(sorted(df_dict.items(), key=lambda item: item[1], reverse=True)[:15])
    top_10_values = {key: int(value) for key, value in top_10_values.items()}
    top_10_values = {key: value for key, value in top_10_values.items() if value != 0}

    plt.figure(figsize=(12, 9))  # Increase figure size for better visualization
    bars = plt.bar(top_10_values.keys(), top_10_values.values(), color='#6baed6')

    plt.grid(axis='y', linestyle='-', alpha=0.5, color='gray')
    plt.gca().set_axisbelow(True)  # Move grid lines behind bars
    plt.xlabel(xaxis_label, fontsize=22, labelpad=15)
    plt.ylabel(yaxis_label, fontsize=22)
    plt.xticks(rotation=45, ha='right', fontsize=20)
    plt.yticks(fontsize=20)

    # Show the box around the plot
    plt.gca().spines['top'].set_visible(True)
    plt.gca().spines['right'].set_visible(True)
    plt.gca().spines['left'].set_linewidth(0.5)  # Adjust line width for left axis
    plt.gca().spines['bottom'].set_linewidth(0.5)  # Adjust line width for bottom axis

    plt.yscale('log') 
    plt.gca().yaxis.set_major_formatter(ScalarFormatter())
    plt.gca().tick_params(axis='y', which='major', pad=15)  # Increase padding for better readability

    plt.tight_layout()  # Adjust layout for better spacing
    plt.savefig(os.path.join(config.base_directory, 'pdb_api_search/analysis/Outputs/Graphs', f'log_{title}.png'), dpi=600)
    plt.show()

def produce_summary_spreadsheet(ptm_ccd_polymer_dict, ptm_ccd_cov_linked_dict, ptm_list, ptm_list_cov_linked, glycan_dict):

    selected_columns = ["ID", "target_residue", "ChEBI", "RESID", "CCD", "PDB_count"]

    ptm_list['PDB_count'] = ptm_list['CCD_combined'].str.split(', ').explode().map(ptm_ccd_polymer_dict).groupby(level=0).sum()
    ptm_list.rename(columns={'CCD_combined': 'CCD'}, inplace=True)
    ptm_list_cov_linked['PDB_count'] = ptm_list_cov_linked['CCD'].str.split(', ').explode().map(ptm_ccd_cov_linked_dict).groupby(level=0).sum()

    new_ptm_list = ptm_list[selected_columns]
    new_ptm_list_cov_linked_list = ptm_list_cov_linked[selected_columns]

    combined_ptm_list = pd.concat([new_ptm_list, new_ptm_list_cov_linked_list], ignore_index=True)  
   
    # Deal with cases where the ID is in both the ptm_list and ptm_cov_linked_list (merge them)
    agg_funcs = {
    'CCD': lambda x: ', '.join(x),
    'PDB_count': 'sum'
    }

    merged_df = combined_ptm_list.groupby('ID').agg(agg_funcs).reset_index()
    combined_ptm_list.set_index('ID', inplace=True)  # Set "ID" as the index for both DataFrames
    merged_df.set_index('ID', inplace=True)
    combined_ptm_list.update(merged_df)
    combined_ptm_list.reset_index(inplace=True)
    combined_ptm_list.drop_duplicates(subset='ID', keep='first', inplace=True)

    glycan_data = pd.DataFrame(glycan_dict.items(), columns=['ID', 'PDB_count'])

    merged_ptm_list = pd.concat([combined_ptm_list, glycan_data], ignore_index=True)
    merged_ptm_list.to_csv(os.path.join(config.base_directory, 'pdb_api_search/analysis/Outputs/PTM_list_with_PDB_count.csv'))

def find_all_pdb_entries(folder_path):
    all_entries = []

    for file_name in os.listdir(folder_path):
        if file_name.endswith('.csv'):
            file_path = os.path.join(folder_path, file_name)
            with open(file_path, 'r', newline='') as csv_file:
                csv_reader = csv.reader(csv_file)
                next(csv_reader) # Skip the header
                for row in csv_reader:
                    all_entries.extend(row)
    
    return all_entries

def get_total_pdb_structures():
    folder_path_cov_linked = os.path.join(config.base_directory, 'pdb_api_search/outputs/ptm_ccd_cov_linked_pdb_search_output')
    folder_path_polymer = os.path.join(config.base_directory, 'pdb_api_search/outputs/ptm_ccd_polymer_pdb_search_output')

    polymer_entries_list = find_all_pdb_entries(folder_path_polymer)
    cov_linked_entries_list = find_all_pdb_entries(folder_path_cov_linked)

    all_pdb_entries_list = polymer_entries_list + cov_linked_entries_list
    all_pdb_entries_set = set(all_pdb_entries_list)

    glycan_data = get_glycan_dict()

    glycan_pdb_entries = sum(glycan_data.values())
    
    total_including_glycans = glycan_pdb_entries + len(all_pdb_entries_set)
    print(total_including_glycans)

if __name__ == '__main__':

    # Dataframes
    keywords_to_add = pd.read_csv(os.path.join(config.base_directory, 'pdb_api_search/inputs/keywords_to_add.csv'), index_col=0)
    
    ptm_list = pd.read_csv(os.path.join(config.base_directory, 'pdb_api_search/outputs/final_ptm_list_w_all_ccd.csv'))
    ptm_list = ptm_list[ptm_list['CCD_combined'].notna()]
    ptm_list = ptm_list.astype(str)
    lipid_polymer_keywords = ['GPI-anchor.', 'Lipoprotein.', 'Myristate.', 'Palmitate.', 'Prenylation.']
    ptm_list['Lipid'] = ptm_list['PTM_type'].apply(lambda x: 'Yes' if any(keyword in x for keyword in lipid_polymer_keywords) else 'No')
    ptm_list.set_index('ID', inplace=True)
    keywords_to_add.set_index('ID', inplace=True)
    ptm_list.update(keywords_to_add)
    ptm_list.reset_index(inplace=True)

    ptm_list_cov_linked = pd.read_csv(os.path.join(config.base_directory, 'pdb_api_search/inputs/cov_linked_ccds.csv'))
    removed_entires_list = pd.read_csv(os.path.join(config.base_directory, 'pdb_api_search/outputs/removed_entires_list_w_ccd.csv'))

    # PDB search output directories
    ptm_ccd_polymer_search_output_directory = os.path.join(config.base_directory, 'pdb_api_search/outputs/ptm_ccd_polymer_pdb_search_output')
    ptm_cov_linked_polymer_search_output_directory = os.path.join(config.base_directory, 'pdb_api_search/outputs/ptm_ccd_cov_linked_pdb_search_output')
    removed_entries_ccd_polymer_search_output_directory = os.path.join(config.base_directory, 'pdb_api_search/outputs/removed_entries_ccd_polymer_pdb_search_output')
    removed_entries_ccd_cov_linked_search_output_directory = os.path.join(config.base_directory, 'pdb_api_search/outputs/removed_entries_cov_linked_polymer_pdb_search_output')
    
    # Outputs dicts
    ptm_ccd_polymer_dict = pdb_value_counts_dictionary(ptm_ccd_polymer_search_output_directory)
    ptm_ccd_cov_linked_dict = pdb_value_counts_dictionary(ptm_cov_linked_polymer_search_output_directory)

    removed_entries_dict = get_removed_entries_dict()
    glycan_dict = get_glycan_dict()
    lipids_dict = get_lipids_dict(ptm_ccd_polymer_dict, ptm_ccd_cov_linked_dict, ptm_list, ptm_list_cov_linked)
    smPTM_dict = get_smPTM_dict(ptm_ccd_polymer_dict, ptm_ccd_cov_linked_dict, ptm_list, ptm_list_cov_linked)
    total_ptm_dict = get_total_ptm_dict()
    # produce_summary_spreadsheet(ptm_ccd_polymer_dict, ptm_ccd_cov_linked_dict, ptm_list, ptm_list_cov_linked, glycan_dict)
    # get_total_pdb_structures()

    # Graphs
    # graph(removed_entries_dict, 'Removed entries - most common CCD in the PDB', 'CCD code', 'Number of PDB structures')
    # log_graph(glycan_dict, 'Number of structures containing glycans in the PDB', 'Glycosylation type', 'Number of PDB structures')
    # log_graph(lipids_dict, 'Number of structures containing lipid PTMs in the PDB', 'Lipidation type', 'Number of PDB structures')
    # graph(smPTM_dict, 'Most common CCD in the PDB', 'CCD code', 'Number of PDB structures')
    # graph(total_ptm_dict, 'Most common PTM type in the PDB', 'PTM type', 'Number of PDB structures')
    