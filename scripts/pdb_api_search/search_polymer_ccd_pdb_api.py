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

def search_ccd_pdb(complete_ccd_file, output_directory_name):
  complete_ccd_file = complete_ccd_file.astype(str)
  complete_ccd_file['CCD_combined'] = complete_ccd_file['CCD_combined'].apply(lambda x: x.split(', '))

  if "CCD_combined" in complete_ccd_file.columns:
      for ccd_list in complete_ccd_file["CCD_combined"]:
          for ccd in ccd_list:
              ccd = str(ccd)
              results = search_pdb_api(ccd)
              if results:
                  result_df = pd.DataFrame({"PDB_ID": results})
                  output_dir = os.path.join(config.base_directory, 'pdb_api_search/outputs', f"{output_directory_name}")

                  os.makedirs(output_dir, exist_ok=True)
                  output_filename = os.path.join(output_dir, f"{ccd}.csv")
                  
                  # Write DataFrame to CSV
                  result_df.to_csv(output_filename, index=False)
    
def search_pdb_api(ccd):
    base_url = "https://search.rcsb.org/rcsbsearch/v2/query"

    search_params = {
      "query": {
        "type": "terminal",
        "label": "text",
        "service": "text",
        "parameters": {
          "attribute": "rcsb_polymer_entity_container_identifiers.chem_comp_monomers",
          "negation": False,
          "operator": "exact_match",
          "value": ccd
        }
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
  ptm_list = pd.read_csv(os.path.join(config.base_directory, 'pdb_api_search/outputs/final_ptm_list_w_all_ccd.csv'))
  removed_entires = pd.read_csv(os.path.join(config.base_directory, 'pdb_api_search/outputs/removed_entires_list_w_ccd.csv'))

  search_ccd_pdb(ptm_list, "ptm_ccd_polymer_pdb_search_output")
  search_ccd_pdb(removed_entires, "removed_entries_ccd_polymer_pdb_search_output")


