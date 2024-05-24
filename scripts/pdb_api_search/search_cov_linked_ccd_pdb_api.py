import requests
import os
import pandas as pd
import sys
import inspect

# Import config.py
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

import config

def search_ccd_pdb(cov_linked_ccd, output_directory_name):
    cov_linked_ccd = cov_linked_ccd.astype(str)
      
    if "CCD" in cov_linked_ccd.columns:
      for ccd in cov_linked_ccd['CCD']:
          ccd = str(ccd)
          results = search_pdb_api(ccd)

          if results:
            result_df = pd.DataFrame({"PDB_ID": results})
            output_dir = os.path.join(config.base_directory, 'pdb_api_search/outputs', f"{output_directory_name}")
            os.makedirs(output_dir, exist_ok=True)
            output_filename = os.path.join(output_dir, f"{ccd}.csv")
            
            result_df.to_csv(output_filename, index=False)
  
def search_pdb_api(ccd):
    base_url = "https://search.rcsb.org/rcsbsearch/v2/query"

    search_params = {
    "query": {
      "type": "group",
      "nodes": [
        {
          "type": "terminal",
          "service": "text",
          "parameters": {
            "attribute": "rcsb_nonpolymer_instance_feature_summary.comp_id",
            "operator": "exact_match",
            "value": ccd
          }
        },
        {
          "type": "terminal",
          "service": "text",
          "parameters": {
            "attribute": "rcsb_nonpolymer_instance_feature_summary.type",
            "operator": "exact_match",
            "value": "HAS_COVALENT_LINKAGE"
          }
        },
        {
          "type": "terminal",
          "service": "text",
          "parameters": {
            "attribute": "rcsb_nonpolymer_instance_feature_summary.count",
            "value": 0,
            "operator": "greater"
          }
        }
      ],
      "logical_operator": "and",
      "label": "text"
    },
            "request_options": {
                "return_all_hits": True
            },
            "return_type": "entry"
    }

    print("Query params:", search_params)

    response = requests.post(base_url, json=search_params)

    if response.status_code == 200: # request was succesful
        data = response.json()
        pdb_ids = [entry["identifier"] for entry in data["result_set"]]
        print("Response:", pdb_ids)

        return pdb_ids
    else:
        print(f"Error: {response.status_code}")
        print("Response:", response.text)

if __name__ == "__main__":
  ptm_cov_linked_ccds = pd.read_csv(os.path.join(config.base_directory, 'pdb_api_search/inputs/cov_linked_ccds.csv'))
  removed_entries_cov_linked_ccds = pd.read_csv(os.path.join(config.base_directory, 'pdb_api_search/inputs/cov_linked_ccds_removed_entries.csv'))

  search_ccd_pdb(ptm_cov_linked_ccds, "ptm_ccd_cov_linked_pdb_search_output")
  search_ccd_pdb(removed_entries_cov_linked_ccds, "removed_entries_cov_linked_polymer_pdb_search_output")
