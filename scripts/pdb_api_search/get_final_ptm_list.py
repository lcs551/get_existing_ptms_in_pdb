import pandas as pd
import os
import sys
import inspect

# Import config.py
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

import config

def get_all_ccds():
        
     # Merge the CCDs identified by resid processing to the CCDs identified with substructure search
    merged_ccd_ptm_list = pd.merge(ptm_list_w_substructure, resid_ccd_match_manual, left_on='RESID', right_on='RESID', how='left')
    merged_ccd_ptm_list['CCD_resid'] = merged_ccd_ptm_list['Check'] 
    merged_ccd_ptm_list.drop(['CCD_ID'], axis=1, inplace=True)
    ptm_list_CCD_combined = merged_ccd_ptm_list.astype(str)
    ptm_list_CCD_combined['CCD_combined'] = ptm_list_CCD_combined.apply(lambda row: ', '.join({x.strip() for col in ['CCD_substructure', 'CCD_resid'] 
                                                                for x in row[col].split(',') if isinstance(row[col], str) and row[col] != 'nan'}), axis=1)
    
    # Merge the manually identified CCDs
    manually_identified_ccd['CCD'] = manually_identified_ccd['CCD'].str.strip()
    ptm_list_CCD_combined['CCD_combined'] = ptm_list_CCD_combined['CCD_combined'].str.strip()

    for index, row in manually_identified_ccd.iterrows():
        # Find the corresponding row in ptm_list_CCD_combined
        idx = ptm_list_CCD_combined.index[ptm_list_CCD_combined['ID'] == row['ID']]
        
        # If the row exists, update CCD_combined
        if len(idx) > 0:
            # Check if the existing value is empty
            if ptm_list_CCD_combined.loc[idx, 'CCD_combined'].iloc[0] == '':
                ptm_list_CCD_combined.loc[idx, 'CCD_combined'] += f"{row['CCD']}"
            else:
                ptm_list_CCD_combined.loc[idx, 'CCD_combined'] += f", {row['CCD']}"
    
    return ptm_list_CCD_combined

if __name__ == "__main__":
        
    resid_ccd_match_manual = pd.read_csv(os.path.join(config.base_directory, 'resid_pdb_search/outputs/extracted_ccd_from_resid_manual_check.csv'))
    ptm_list_w_substructure = pd.read_csv(os.path.join(config.base_directory, 'substructure_search/outputs/ptm_list_w_substructure_ccd.csv'))
    manually_identified_ccd = pd.read_csv(os.path.join(config.base_directory, 'pdb_api_search/inputs/manually_identified_ccds.csv'))
    include_ptm_list = pd.read_csv(os.path.join(config.base_directory, 'pdb_api_search/inputs/processed_uniprot_ptm_list_include.csv'))

    # Get all CCDs, then separate into PTM list and removed entries (non-PTMs, such as chromophores, metaboiltes etc.)
    all_ccd_ptm_list = get_all_ccds()

    yes_ptm_list = include_ptm_list[include_ptm_list['Include'] == 'Yes']['ID']
    no_ptm_list = include_ptm_list[include_ptm_list['Include'] == 'No']['ID']

    final_ptm_list = all_ccd_ptm_list[all_ccd_ptm_list['ID'].isin(yes_ptm_list)]
    removed_ptm_list = all_ccd_ptm_list[all_ccd_ptm_list['ID'].isin(no_ptm_list)]

    final_ptm_list.to_csv(os.path.join(config.base_directory, 'pdb_api_search/outputs/final_ptm_list_w_all_ccd.csv'), index=False)
    removed_ptm_list.to_csv(os.path.join(config.base_directory, 'pdb_api_search/outputs/removed_entires_list_w_ccd.csv'), index=False)
