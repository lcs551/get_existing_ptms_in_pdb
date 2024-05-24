from gemmi import cif
import os
import pandas as pd
from tqdm import tqdm
import sys
import inspect

# Import config.py
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

import config

def get_mod_res_covalently_bound(doc):
    block = doc.sole_block()

    ### search struct_conn.ptnr1 and 2
    # found_mod_res = False
      
    partner_1 = []
    partner_2 = []
    connection_type = []

    column = block.find_loop('_struct_conn.ptnr1_label_comp_id')
    for value in column:
        partner_1.append(value)
    
    column = block.find_loop('_struct_conn.ptnr2_label_comp_id')
    for value in column:
        partner_2.append(value)

    column = block.find_loop('_struct_conn.conn_type_id')
    for value in column:
        connection_type.append(value)

    partner_1_bond_type = dict(zip(partner_1, connection_type))
    partner_2_bond_type = dict(zip(partner_2, connection_type))
    
    partner_1_covale_dict = {key: value for key, value in partner_1_bond_type.items() if value == 'covale'}
    partner_2_covale_dict = {key: value for key, value in partner_2_bond_type.items() if value == 'covale'}
    
    covalent_partners = {**partner_1_covale_dict, **partner_2_covale_dict}

     ### get modified residue ###
    mod_res = []

    column = block.find_loop('_pdbx_struct_mod_residue.label_comp_id')
    for key, value in covalent_partners.items():
        if value in covalent_partners:
            for val in column:
                mod_res.append(val)
    mod_res = list(set(column)) # only get unique values
  
    return mod_res

def read_csv_files(directory):
    data = {}
    for filename in os.listdir(directory):
        if filename.endswith(".csv"):
            file_path = os.path.join(directory, filename)
            df = pd.read_csv(file_path, header=None, skiprows=1)
            key = os.path.splitext(filename)[0]
            data[key] = df[0].tolist()
    return data

if __name__ == "__main__":
    directory = "/vault/pdb_mirror/data/structures/all/mmCIF/" # Location of PDB mirror (CHANGE)
    directory_path =  os.path.join(config.base_directory, 'resid_pdb_search/outputs/resid_pdb_search_outputs')
    resid_dict = {}
    
    result_dict = read_csv_files(directory_path)
    # total_length = sum(len(value) for value in result_dict.values())

    for key,value in tqdm(result_dict.items()):
        ccd_dict = {}
        for pdb in value:
            file_path = os.path.join(directory, f"{pdb.lower()}.cif.gz")
            doc = cif.read(file_path)
            ccd = get_mod_res_covalently_bound(doc)

            if ccd:
                file_name = os.path.splitext(pdb)[0]
                ccd_dict[file_name] = ccd
        
        resid_dict[key]=ccd_dict

        output_path = os.path.join(config.base_directory, 'resid_pdb_search/outputs/resid_ccd_outputs', f'{key}.csv')
        with open(output_path, 'w') as output_file:
            output_file.write("PDB_ID,CCD_ID\n")
            for pdb,ccds in ccd_dict.items():
                output_file.write(f"{pdb},{'/'.join(ccds)}\n")
