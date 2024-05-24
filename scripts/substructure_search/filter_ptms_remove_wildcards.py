from rdkit import Chem
import pandas as pd
import sys
import os
import inspect

# Import config.py
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

import config

def filter_sdf_by_chebi(sdf_file, chebi_ids, output_file):
    writer = Chem.SDWriter(output_file)
    suppl = Chem.SDMolSupplier(sdf_file)

    for mol in suppl:
        if mol is None:
            continue
        props = mol.GetPropsAsDict()
        chebi_id = props.get('ChEBI ID')  # ChEBI ID is stored as a property

        if chebi_id in chebi_ids:
            writer.write(mol)
    writer.close()

def get_ptm_chebi_sdf():
    uniprot_ptm_file = pd.read_csv(os.path.join(config.base_directory, 'uniprot_ptm_list/outputs/final_uniprot_ptm_list.csv'), dtype=str)
    uniprot_ptm_file = uniprot_ptm_file.dropna(subset=["ChEBI"])
    uniprot_ptm_file['ChEBI'] = 'CHEBI:' + uniprot_ptm_file['ChEBI'].astype(str)
    chebi_ids_list = uniprot_ptm_file['ChEBI'].tolist()
    chebi_ids_list = [chebi_id.rstrip('.0') for chebi_id in chebi_ids_list]
    
    chebi_sdf = os.path.join(config.chebi_sdf_file_path)
    output_file = os.path.join(config.base_directory, 'substructure_search/outputs/chebi_ptms.sdf')
    filter_sdf_by_chebi(chebi_sdf, chebi_ids_list, output_file)

# Set carbon's wildcard atom to O, and remove wildcard atom and bond from N
def modify_molecule(mol):
    try:
        # Aromatize the molecule
        Chem.SanitizeMol(mol, Chem.SANITIZE_ALL, catchErrors=True)
        Chem.Kekulize(mol, clearAromaticFlags=True)

        # Find wildcard atoms bonded to carbon
        for bond in mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            if mol.GetAtomWithIdx(begin_idx).GetSymbol() == '*' and mol.GetAtomWithIdx(end_idx).GetSymbol() == 'C':
                mol.GetAtomWithIdx(begin_idx).SetAtomicNum(8)  # Set the wildcard atom to Oxygen
                mol.GetAtomWithIdx(begin_idx).SetFormalCharge(0)  # Set formal charge to 0
            elif mol.GetAtomWithIdx(end_idx).GetSymbol() == '*' and mol.GetAtomWithIdx(begin_idx).GetSymbol() == 'C':
                mol.GetAtomWithIdx(end_idx).SetAtomicNum(8)  # Set the wildcard atom to Oxygen
                mol.GetAtomWithIdx(end_idx).SetFormalCharge(0)  # Set formal charge to 0

        # Create a mutable version of the molecule to modify bonds and atoms directly
        mol = Chem.RWMol(mol)

        # Remove wildcard atoms bonded to N atoms
        to_remove = []
        for bond in mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            if (mol.GetAtomWithIdx(begin_idx).GetSymbol() == '*' and
                    mol.GetAtomWithIdx(end_idx).GetSymbol() == 'N') or \
                    (mol.GetAtomWithIdx(begin_idx).GetSymbol() == 'N' and
                     mol.GetAtomWithIdx(end_idx).GetSymbol() == '*'):
                to_remove.append((begin_idx, end_idx))

        # Remove bonds and wildcard atoms
        for begin_idx, end_idx in to_remove:
            mol.RemoveBond(begin_idx, end_idx)
            if mol.GetAtomWithIdx(begin_idx).GetSymbol() == '*':
                mol.RemoveAtom(begin_idx)
            elif mol.GetAtomWithIdx(end_idx).GetSymbol() == '*':
                mol.RemoveAtom(end_idx)

        for atom in mol.GetAtoms():
            atom.SetFormalCharge(0)

        return mol
    except Exception as e:
        print(f"Error in modifying molecule: {e}")
        return None

if __name__ == "__main__":
    
    get_ptm_chebi_sdf()

    input_sdf_file = os.path.join(config.base_directory, 'substructure_search/outputs/chebi_ptms.sdf')
    output_sdf_file = os.path.join(config.base_directory, 'substructure_search/outputs/final_chebi_ptms.sdf')

    suppl = Chem.SDMolSupplier(input_sdf_file)
    writer = Chem.SDWriter(output_sdf_file)
    num_structures = 0

    # Iterate over each molecule in the input file, modify it, and write to the output file
    for idx, mol in enumerate(suppl):
        if mol is not None:
            modified_mol = modify_molecule(mol)
            if modified_mol is not None:
                writer.write(modified_mol)

    writer.close()
    print("Conversion complete.")
