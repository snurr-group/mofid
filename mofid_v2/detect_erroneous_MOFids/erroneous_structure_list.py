import pandas as pd
from pathlib import Path
import glob

CORE_PATH = ".../core2024-v5"

def get_folder_dir(cifname):

    """Assumes a certain directory structure to fetch the directory of the given filename."""
    dir_p = {'ASR':"dataset-ASR", 'FSR':"dataset-FSR", 'ion':"dataset-ion"}
    cifname = cifname.split('_pacman')[0]+"_pacman"
    type = cifname.split('_pacman')[0].split('_')[-1]
    import re
    cif_pattern = re.compile(r'^[a-zA-Z].*$')
    if cif_pattern.match(cifname):
        subfolder_name = cifname[0]
    else:
        subfolder_name = 'others'
    base_dir = Path(f"{CORE_PATH}/dataset-{type}/{subfolder_name}/Output/{cifname}/")
    if not base_dir.is_dir():
        print(f"WARNING - base_dir is not a dir! {base_dir}")

    # print(f"debug:{base_dir}")
    return base_dir   # This might need to be adjusted depending on the actual structure


def read_and_combine_datasets():
    # Read the CSV files with multi-level columns
    df_list = []
    for name in ["FSR", "ASR", "ion"]:
        df = pd.read_csv(f"../result_data/dataset-{name}_solvent_check.csv", header=[0, 1])
        df_list.append(df)
    combined_df = pd.concat(df_list, ignore_index=True)
    combined_df = combined_df.replace('Yes [1]', 'Yes')
    # breakpoint()
    
    # Read and merge MOFID datasets
    df_mofid_list = []
    for name in ["FSR", "ASR", "ion"]:
        df = pd.read_csv(f"../result_data/dataset-{name}_mofid.csv")
        df_mofid_list.append(df)
    combined_df_mofid = pd.concat(df_mofid_list, ignore_index=True)
    
    return combined_df, combined_df_mofid

def filter_type1_structures(df):
    # Filtering logic for Type 1 structures
    condition = (df[("bound_solvent", "MetalOxo")] != "No") | (df[("free_solvent", "MetalOxo")] != "No")

    # Filter the DataFrame based on the condition
    filtered_df = df[condition]
    
    # Select the specific columns needed along with filtering
    # Assuming 'mof name' is the label of the first column in a non-multiindex way
    result_df = filtered_df.loc[:, [("Solvent_Type","Algorithm"), ("bound_solvent", "MetalOxo"), ("free_solvent", "MetalOxo")]].copy()
    
    # Ensure 'mof name' is the first column if it's not already a DataFrame default
    # If 'mof name' needs to be corrected or specifically set, handle that here
    result_df.columns = ['mof_name', 'bound_solvent_MetalOxo', 'free_solvent_MetalOxo']

    return result_df

def determine_type(mof_name):
    # Determine the type based on mof_name
    if "_ASR_" in mof_name:
        return "ASR"
    elif "_FSR_" in mof_name:
        return "FSR"
    elif "_ion_" in mof_name:
        return "ion"
    return "Unknown"

def get_files_as_string(mof_name, pattern):
    # Use glob to search for files matching pattern
    folder_dir = get_folder_dir(mof_name)
    file_list = [x.name for x in list(folder_dir.glob(f"{mof_name}_{pattern}_*.xyz"))]
    return str(file_list)

def main():
    solvent_check, mofid = read_and_combine_datasets()

    # Type 1 DataFrame
    type1_df = filter_type1_structures(solvent_check)
    type1_df['mofid'] = type1_df['mof_name'].map(mofid.set_index('cifname')['mofid'])
    type1_df['type'] = type1_df['mof_name'].apply(determine_type)

    # Type 2 DataFrame logic here, similar to Type 1 but with different conditions
    condition2a = (solvent_check[("bound_solvent", "MetalOxo")] == "No") & (solvent_check[("bound_solvent", "AllNode")] == "Yes")
    condition2b = (solvent_check[("free_solvent", "MetalOxo")] == "No") & (solvent_check[("free_solvent", "AllNode")] == "Yes")
    type2_df = solvent_check[condition2a | condition2b]
    type2_df = type2_df.loc[:, [("Solvent_Type","Algorithm")]].copy()
    type2_df.columns = ['mof_name']
    type2_df['bound_solvent_xyz'] = type2_df['mof_name'].apply(lambda x: get_files_as_string(x, "AllNode_bound"))
    type2_df['free_solvent_xyz'] = type2_df['mof_name'].apply(lambda x: get_files_as_string(x, "AllNode_free"))
    type2_df['mofid'] = type2_df['mof_name'].map(mofid.set_index('cifname')['mofid'])
    type2_df['type'] = type2_df['mof_name'].apply(determine_type)

    type_order = ["ASR", "FSR", "ion"]

    # Convert the 'type' column to a categorical type with a defined order
    type2_df['type'] = pd.Categorical(type2_df['type'], categories=type_order, ordered=True)
    type1_df['type'] = pd.Categorical(type1_df['type'], categories=type_order, ordered=True)

    # Sort the DataFrame by the 'type' column
    type1_df = type1_df.sort_values(by='type')
    type2_df = type2_df.sort_values(by='type')

    # Save to CSV or handle further as needed
    type1_df.to_csv('Type_1_Structures.csv', index=False)
    type2_df.to_csv('Type_2_Structures.csv', index=False)

if __name__ == "__main__":
    main()
