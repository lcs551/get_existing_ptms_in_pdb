import requests
import pandas as pd
import os
import sys
import inspect

# Import config.py
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

import config

def search_pdb_api(resid):
    base_url = "https://search.rcsb.org/rcsbsearch/v2/query"

    search_params = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text_chem",
                    "parameters": {
                        "attribute": "rcsb_chem_comp_related.resource_accession_code",
                        "operator": "exact_match",
                        "negation": False,
                        "value": resid
                    }
                },
                {
                    "type": "terminal",
                    "service": "text_chem",
                    "parameters": {
                        "attribute": "rcsb_chem_comp_related.resource_name",
                        "operator": "exact_match",
                        "value": "RESID",
                        "negation": False
                    }
                }
            ],
            "label": "nested-attribute"
        },
        "request_options": {
            "return_all_hits": True
        },
        "return_type": "entry"
    }

    print("Query params:", search_params)  # Print query params for debugging

    response = requests.post(base_url, json=search_params)

    # Check if the request was successful (status code 200)
    if response.status_code == 200:
        # Parse the JSON response and extract PDB IDs
        data = response.json()
        pdb_ids = [entry["identifier"] for entry in data["result_set"]]
        
        print("Response:", pdb_ids)  # Print response for debugging

        return pdb_ids
    else:
        print(f"Error: {response.status_code}")
        print("Response:", response.text)  # Print response for debugging

if __name__ == "__main__":
    uniprot_ptm_file = pd.read_csv(os.path.join(config.base_directory, 'uniprot_ptm_list/outputs/final_uniprot_ptm_list.csv'))
    standard_amino_acid_resid_list = pd.read_csv(os.path.join(config.base_directory, 'resid_pdb_search/inputs/resid_saa.csv'))

    if "RESID" in uniprot_ptm_file.columns:
        for resid in uniprot_ptm_file["RESID"]:
            if not pd.isnull(resid):
                resid = str(resid)
                if resid not in standard_amino_acid_resid_list['RESID'].values:
                    results = search_pdb_api(resid)
                    if results:  # Check if results list is not empty
                        # Create DataFrame with a column named 'PDB_ID'
                        result_df = pd.DataFrame({"PDB_ID": results})
                        
                        output_dir = os.path.join(config.base_directory, 'resid_pdb_search/outputs/resid_pdb_search_outputs')
                        os.makedirs(output_dir, exist_ok=True)
                        output_filename = os.path.join(output_dir, f"{resid}.csv")
                        
                        # Write DataFrame to CSV
                        result_df.to_csv(output_filename, index=False)
            else:
                continue
