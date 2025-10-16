# MOFid-v2 Conversion and Analysis Scripts

## Overview
This repository contains scripts for processing and analyzing CoRE MOF 2024 structures using the MOFid software. The scripts enable the decomposition of CIF files into building blocks, such as nodes and linkers, followed by grouping and further characterization based on their chemical compositions.

## Instructions for Using the Scripts

### Step 1: Extracting Building Blocks
Using the MOFid software, we processed all CIF files from the CoRE MOF 2024 dataset. This process generated:
- **MOFid-v1 identifiers** for the MOFs.
- Decomposition of CIF files into their constituent building blocks (nodes and linkers).

### Step 2: Processing Nodes
1. **Splitting Nodes into Individual `.xyz` Files**:
   - The **AllNode method** from MOFid provides the coordinates of all nodes in a MOF.
   - To simplify analysis, these node files were split into individual `.xyz` files based on their chemical compositions using the script:
     ```bash
     python split_nodes_to_xyz.py
     ```

2. **Grouping Nodes by Composition**:
   - After generating `.xyz` files for all nodes, the nodes were grouped by their chemical compositions using:
     ```bash
     python group_node_by_composition.py
     ```

3. **Characterizing Node Types**:
   - Further analysis was performed to classify nodes into different types based on their structures using:
      ```bash
     python group_node_by_structure.py
     ```  

### Step 3: Processing Linkers
1. **Splitting Linkers into Individual `.xyz` Files**:
   - The **MetalOxo method** from MOFid provides the coordinates of linkers.
   - These linkers were split into individual `.xyz` files based on their chemical compositions using:
     ```bash
     python split_linkers_to_xyz.py
     ```

2. **Grouping Linkers by SMILES String**:
   - Using the output of MOFid, the linkers were grouped by their canonical SMILES representation as implemented in OpenBabel:
     ```bash
     python group_linkers_by_SMILES.py
     ```

3. **Converting Linkers to XYZ Files**:
   - Using OpenBabel, the SMILES strings were reconverted to XYZ format using:
     ```bash
     python convert_smiles_to_xyz.py
     ```

## Notes
- Ensure that all necessary dependencies and input files are correctly configured before running the scripts.
- Refer to the CoRE MOF 2024 publication for detailed descriptions of the dataset and methods.

