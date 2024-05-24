from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import csv
import os
import sys
import inspect

# Import config.py
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

import config

def add_ccd_substructure(input_file, output_file, matches):
    # Create a dictionary to store ChEBI_ID to list of CCD mappings
    chebi_to_ccd = {}
    for ccd_name, chebi_name in matches:
        if chebi_name not in chebi_to_ccd:
            chebi_to_ccd[chebi_name] = []
        chebi_to_ccd[chebi_name].append(ccd_name)

    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames + ['CCD_substructure']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for row in reader:
            chebi_id = row['ChEBI']
            # Add "CHEBI:" prefix if missing from input ChEBI_ID
            chebi_id = "CHEBI:" + chebi_id
            chebi_id = chebi_id.rstrip(".0")
            ccd_substructures = chebi_to_ccd.get(chebi_id, [])  # Get corresponding CCD values or empty list if not found
            if not ccd_substructures:
                ccd_substructures = ['N/A']  # If no CCD found, set it as 'N/A'
            ccd_substructure_str = ', '.join(ccd_substructures)
            row['CCD_substructure'] = ccd_substructure_str
            writer.writerow(row)

if __name__ == "__main__":
    ccd_sdf = Chem.SDMolSupplier(os.path.join(config.base_directory, 'substructure_search/inputs/Components-pub.sdf'))
    ccd_molecules = [mol for mol in ccd_sdf if mol is not None]

    chebi_sdf = Chem.SDMolSupplier(os.path.join(config.ccd_sdf_file_path))    
    chebi_molecules = [mol for mol in chebi_sdf if mol is not None]

    # Substructure search
    matches = []
    for chebi_mol in chebi_molecules:
        for ccd_mol in ccd_molecules:
            # Check if both molecules have the same number of atoms and bonds
            if ccd_mol.GetNumAtoms() == chebi_mol.GetNumAtoms() and ccd_mol.GetNumBonds() == chebi_mol.GetNumBonds():
                match_atoms = ccd_mol.GetSubstructMatch(chebi_mol)
                if match_atoms:
                    print(match_atoms)
                if match_atoms and len(match_atoms) == chebi_mol.GetNumAtoms():
                    matches.append((ccd_mol.GetProp('_Name'), chebi_mol.GetProp('ChEBI ID')))

    input_file = os.path.join(config.base_directory, 'uniprot_ptm_list/outputs/final_uniprot_ptm_list.csv')
    output_file = os.path.join(config.base_directory, 'substructure_search/outputs/ptm_list_w_substructure_ccd.csv')

    add_ccd_substructure(input_file, output_file, matches)

    print("Data written to .csv")
