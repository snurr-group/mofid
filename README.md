# MOFid
A system for rapid identification and analysis of metal-organic frameworks

**WARNING: code is under active development and refactoring. Details subject to change**

## Objective
Supplement the current MOF naming conventions with a canonical, machine-readable identifier to facilitate data mining and searches. Accomplish this goal by representing MOFs according to their nodes + linkers + topology

## Requirements
1. A Python environment is required. If you do not have Python environment installed, we recommend downloading and installing [Anaconda](https://www.anaconda.com/distribution/#download-section). MOFid is compatible with both Py2 and Py3.
2. Make sure you have the following: a C++ compiler, [`cmake`](https://cmake.org/), and access to GNU commands (such as `make`). These are all typically available on Linux machines. If running on Windows, we recommend installing [Cygwin](https://www.cygwin.com/) and including both the `make` and `wget` packages in addition to the default options during the setup process.
3. Make sure you have the [Java Runtime Environment](https://www.java.com/en/download/) installed and included in your system's path, which is needed for [Systre](http://gavrog.org/).

## Installation
1. Run `make init` in the base `/mofid` directory.
2. If using Python 2.x, you'll need to install the [`subprocess32`](https://pypi.org/project/subprocess32/) package, which can be done via `pip install --user subprocess32`.
3. Run `make test` to verify that everything is set up correctly. Some warning messages are expected and normal, because the test suite includes some non-MOF structures. At the end, the test code will report the number of mismatches, which should be zero (except for a known issue with ZIFs on some computers).

## Usage
1. Run `/mofid/Python/run_mofid.py NAME_OF_CIF.cif` to output the unique identifier for one MOF of interest. Decomposition information is written to `/mofid/Output/` by default.
2. Run `bash /mofid/Scripts/run_folder.sh CIF_DIR OUTPUT_DIR > results.smi 2> log.txt` to output the unique identifiers for a folder of MOFs, using a temporary `OUTPUT_DIR` for scratch work.  The list of MOFids and MOFkeys will also be written to files like `OUTPUT_DIR/folder_mofkey.tsv`.
3. Run `/mofid/Python/convert_smi_to_tables.py results.smi` to translate a list of MOFid's from a .smi file into .tsv tables describing the generated SMILES and topological data about the MOFs to the `/mofid/TableOutput/` folder.
