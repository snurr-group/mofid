# MOFid
A system for rapid identification and analysis of metal-organic frameworks

**WARNING: code is under active development and refactoring. Details subject to change**

## Objective
Supplement the current MOF naming conventions with a canonical, machine-readable identifier to facilitate data mining and searches. Accomplish this goal by representing MOFs according to their nodes + linkers + topology

## Requirements
1. A Python environment is required. If you do not have Python environment installed, we recommend downloading and installing [Anaconda](https://www.anaconda.com/distribution/#download-section). MOFid is compatible with both Python 2/3.
2. Make sure you have the following: a C++ compiler, [`cmake`](https://cmake.org/), and access to GNU commands (such as `make`). These are all typically available on Linux machines. If running on Windows, we recommend using [Cygwin](https://www.cygwin.com/) and including both the `make` and `wget` packages in addition to the default options during the Cygwin installation process.
3. Make sure you have the [Java Runtime Environment](https://www.java.com/en/download/) installed and included in your system's path. If unsure, try running `java` in the command line to see if it successfully calls Java.

## Installation
1. Run `make init` in the base `/mofid` directory.
2. Run `python set_paths.py; pip install .` in the base `/mofid` directory.  If you encounter permissions errors (typically not with anaconda), you may need to run `pip install --user .`

## Usage
In a Python script, the user simply has to call the `run_mofid.cif2mofid(cif_path,output_path='Output')` function. The first argument is required and is the path to the MOF CIF. The second argument is optional and is the directory to store the MOFid decomposition information, which defaults to `/Output` if not specified. An example of how to call MOFid is shown below.
```python
from mofid.run_mofid import cif2mofid
cif_path = '/path/to/my/mof.cif'
mofid = cif2mofid(cif_path)
```
The output of the `mofid.cif2mofid` function is a dictionary containing six entries: the MOFid (`mofid`), MOFkey (`mofkey`), SMILES string (`smiles`), topology (`topology`), catenation (`cat`), and basename of the CIF (`cifname`).