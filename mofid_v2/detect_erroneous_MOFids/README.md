## Overview

This repository contains scripts for identifying erroneous cases in the CoRE MOF 2024 structures using the MOFid software. These scripts check for the presence of spurious solvent fragments in the outputs of the MetalOxo and AllNode algorithms for each MOF structure.

## Using the Scripts

The script `analyze_solvent_containing_mofs.py` generates a separate table containing solvent fragment information for each of the three datasets: ASR, FSR, and ion datasets.

The script `erroneous_structure_list.py` generates two lists of MOFs whose MOFid may be incorrect. It uses the MOFid output folders and split fragments to identify the chemical information of solvent fragments. Specify the root path to all three dataset folders using the variable `CORE_PATH`.

The output files are:

- **Type_1_Structures.csv**: This table stores information about MOFs that contain linker parsing issues.
- **Type_2_Structures.csv**: This table stores information about MOFs that contain topology simplification issues.

## Notes

- Ensure that all necessary dependencies and input files are correctly configured before running the scripts.
- Refer to the CoRE MOF 2024 publication for detailed descriptions of the dataset and methods.