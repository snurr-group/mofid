"""
returns a table that documents which folders (corresponding to cif names) contain non-empty solvent files
@author: Kunhuan Liu
"""

from pathlib import Path
import glob
import pandas as pd
from pymatgen.analysis.structure_matcher import StructureMatcher


def check_solvent_folder(target_F):
    """
    return a dictionary of {file_id, data}, where data is 3x2 :
    {MetalOxo, SingleNode, AllNode: {'free_solvent', 'bound_solvent': No/Yes*/Yes**/Yes***}}
    """
    directory_path = Path(target_F)
    dirs = {'MetalOxo': directory_path / 'MetalOxo',
            'SingleNode': directory_path / 'SingleNode',
            'AllNode': directory_path / 'AllNode'}

    summ_dict = {algo: {} for algo in dirs}
    free_structure_bank = []
    bound_structure_bank = []
    matcher = StructureMatcher(scale=False)

    def fill_dict(directory_path,solv, algo, summ_dict, structure_bank):
        solvent_path = directory_path / (solv + ".cif")
        print(f"checking: {solvent_path}")
        check = solvent_file_check(solvent_path)
        if check is None:
            summ_dict[algo][solv] = 'NaN'
        elif check is False:
            summ_dict[algo][solv] = 'No'
        else:
            flag = False
            for i, s in enumerate(structure_bank):
                if matcher.fit(check, s):
                    flag = i
            if flag is False:
                # flag = '*' * (len(structure_bank)+1)
                flag = len(structure_bank)+1
            else:
                # flag = '*' * (i + 1)
                flag = i + 1
            summ_dict[algo][solv] = f'Yes [{flag}]'

    for algo, directory_path in dirs.items():
        fill_dict(directory_path, "free_solvent", algo, summ_dict, free_structure_bank)
        fill_dict(directory_path, "bound_solvent", algo, summ_dict, bound_structure_bank)


    return summ_dict
    

def check_solvent_files(target_F, algorithm="SingleNode", include_bound = False):
    '''
    output: a df with 6 dimensions
    '''


    # Construct the directory path
    directory_path = Path(target_F) / algorithm

    # Check for free_solvent.cif
    free_solvent_path = directory_path / "free_solvent.cif"
    free_solvent_status = solvent_file_check(free_solvent_path)

    
    if include_bound: # return it contains solvent (TRUE) when either has solvent
        bpath = directory_path / "bound_solvent.cif"
        bstatus = solvent_file_check(bpath)
        if free_solvent_status is None and bstatus is None:
            return None
        else:
            return (free_solvent_status or bstatus)
    else:
        return free_solvent_status

def solvent_file_check(file_path):
    try:
        with open(file_path, 'r') as file:
            for i, _ in enumerate(file, 1):
                if i > 21:
                    # contains atoms
                    from pymatgen.io.cif import CifParser
                    parser = CifParser(file_path, site_tolerance = 1e-6, occupancy_tolerance = 1.0)
                    structure = parser.parse_structures()[0]
                    return structure
                    # return True  # File has more than 21 lines, hence contains atoms
            return False  # File has 21 or fewer lines, hence "empty"
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        return None





def get_mofid(target_F):
    try:
        with open(target_F / "python_mofid.txt") as file:
            for line in file:
                if line.strip():
                    return line # line is not empty
    except Exception as e:
        print(f"Encountered an issue when grabbing python_mofid.txt:  {target_F}\n {e}")
    return None




def process_OutputFolder(dest_F, include_bound):
    # Initialize an empty DataFrame
    # df = pd.DataFrame(columns=['file_name', 'SN_solvent', 'AN_solvent'])
    df_data = []
    current_folder = Path(dest_F).name
    output_path = Path(dest_F) / "Output"
    folders = sorted([f for f in output_path.glob('*') if f.is_dir()])

    for folder in folders:
        folder_path = Path(folder)
        file_name = folder_path.name

        # Check solvents for SingleNode
        sn_result = check_solvent_files(folder, "SingleNode", include_bound)
        sn_status = interpret_result(sn_result)

        # Check solvents for AllNode
        an_result = check_solvent_files(folder, "AllNode", include_bound)
        an_status = interpret_result(an_result)

        # Grab MOFid
        mofid = get_mofid(folder)

        # Append to DataFrame
        df_data.append({'file_name': file_name, 'SN_solvent': sn_status, 'AN_solvent': an_status, 'mofid' : mofid, 'folder': current_folder})

    # df = pd.DataFrame(df_data)
    return df_data


def masking_mofid(original_mofid,reason = "linker_fragmentation_error"):
    # decompose original mofid
    filename = original_mofid.split(';')[1:-1]
    return f"* MOFid-v2.NA.NA;{filename};{reason}"



import pandas as pd
from pathlib import Path

def process_structure_folders(folder_list):
    # Initialize the list that will hold all data rows
    data = []

    # Iterate through each folder in the provided list
    for folder_name in folder_list:
        # Call your existing function
        result_dict = check_solvent_folder(folder_name)
        
        # Flatten the dictionary structure to fit the DataFrame format with simplified column names
        row_data = {
            f"{algo}_{suffix[0]}": value
            for algo, subdict in result_dict.items()
            for suffix, value in subdict.items()
        }
        
        # Add folder name or structure name to the row_data
        row_data['Structure'] = Path(folder_name).name
        
        # Append the row to the data list
        data.append(row_data)

    # Convert the list of dictionaries to a DataFrame
    df = pd.DataFrame(data)

    # Reorder the DataFrame columns to start with 'Structure' followed by the data columns
    first_col = ['Structure']
    other_cols = sorted([col for col in df.columns if col != 'Structure'])
    df = df[first_col + other_cols]

    return df

def process_structure_folders_multiindex(folder_list):
    # Initialize the list that will hold all data rows
    data = []

    # Iterate through each folder in the provided list
    for folder in folder_list:
        # Call your existing function
        result_dict = check_solvent_folder(folder)
        folder_name = Path(folder).name
        
        # Prepare row as a flat list with multi-indexed columns in mind
        row_data = []
        column_index = []
        
        for algo, subdict in result_dict.items():
            for solvent_type, value in subdict.items():
                # Append the value to row data
                row_data.append(value)
                # Append a tuple (algorithm, solvent_type) to the column index list
                column_index.append((algo, solvent_type))
        
        # Append a tuple of the row data and the folder name to the data list
        data.append((folder_name, row_data))

    # Create a DataFrame using a MultiIndex for the columns
    # multi_index = pd.MultiIndex.from_tuples(column_index, names=['Solvent_Type', 'Algorithm'])
    multi_index = pd.MultiIndex.from_tuples(column_index, names=['Algorithm', 'Solvent_Type'])
    df = pd.DataFrame([row_data for _, row_data in data], columns=multi_index, index=[folder_name for folder_name, _ in data])

    df = df.swaplevel('Solvent_Type', 'Algorithm', axis=1)
    df.sort_index(axis=1, inplace=True)
    algorithm_order = ['MetalOxo', 'SingleNode', 'AllNode']
    desired_index = pd.Index(algorithm_order, name='Algorithm')
    df_reordered = df.reindex(columns=desired_index, level='Algorithm')

    return df_reordered


def main():
    root_dirs = ['dataset-ion', 'dataset-ASR', 'dataset-FSR']

    for root_dir in root_dirs:
        # Create a Path object for the root directory
        path = Path(root_dir)
        
        # Grab all folders that match '*/Output/*/'
        subfolder_list = [folder for folder in path.glob('*/Output/*/') if folder.is_dir()]

        # Call the function to process these folders and get a DataFrame
        df = process_structure_folders_multiindex(subfolder_list)

        # Save the DataFrame to a CSV file
        df.to_csv(f"{root_dir}_solvent_check.csv")

if __name__ == "__main__":

    main()
