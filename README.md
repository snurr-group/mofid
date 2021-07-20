# MOFid
A system for rapid identification and analysis of metal-organic frameworks.

Check out [https://snurr-group.github.io/web-mofid/](https://snurr-group.github.io/web-mofid/) to quickly and easily run these tools in your browser! No programming skills required. If you wish to generate the MOFid for a large number of structures, see the Python-based interface described below.

Please cite [DOI: 10.1021/acs.cgd.9b01050](https://pubs.acs.org/doi/abs/10.1021/acs.cgd.9b01050) if you use MOFid in your work.

## Objective
Supplement the current MOF naming conventions with a canonical, machine-readable identifier to facilitate data mining and searches. Accomplish this goal by representing MOFs according to their nodes + linkers + topology

## Requirements
1. A Python environment is required. If you do not have a Python environment installed, we recommend downloading and installing [Anaconda](https://www.anaconda.com/distribution/#download-section). MOFid is compatible with both Python 2/3.
2. Make sure you have the following: a C++ compiler, [`cmake`](https://cmake.org/), and access to GNU commands (such as `make`). These are all typically available on Linux machines. If running on Windows, we recommend using [Cygwin](https://www.cygwin.com/) and including the `cmake`, `make`, `wget`, `gcc-core`, `gcc-g++`, and `pkg-config` packages in addition to the default options during the Cygwin installation process.
3. Make sure you have the [Java Runtime Environment](https://www.java.com/en/download/) installed and included in your system's path. If unsure, try running `java` in the command line to see if it successfully calls Java.

## Installation
1. Run `make init` in the base `/mofid` directory.
2. Run `python set_paths.py` followed by `pip install .` in the base `/mofid` directory.  If you encounter permissions errors (typically not with Anaconda), you may need to run `pip install --user .`

## Usage
In a Python script, the user simply has to call the `run_mofid.cif2mofid(cif_path,output_path='Output')` function. The first argument is required and is the path to the MOF CIF. The second argument is optional and is the directory to store the MOFid decomposition information, which defaults to `/Output` if not specified. An example of how to call MOFid is shown below.
```python
from mofid.run_mofid import cif2mofid
cif_path = '/path/to/my/mof.cif'
mofid = cif2mofid(cif_path)
```
The output of the `mofid.cif2mofid` function is a dictionary containing eight entries: the MOFid (`mofid`), MOFkey (`mofkey`), SMILES string (`smiles`, `smiles_nodes`, or `smiles_linkers`), topology (`topology`), catenation (`cat`), and basename of the CIF (`cifname`).

## Background and Troubleshooting
Please read the page [here](https://github.com/snurr-group/web-mofid/blob/master/README.md) for a detailed background and for important tips/tricks to help troubleshoot any problematic scenarios.

## Credits
This work is supported by the U.S. Department of Energy, Office of Basic 
Energy Sciences, Division of Chemical Sciences, Geosciences and 
Biosciences through the Nanoporous Materials Genome Center under award 
DE-FG02-17ER16362.

The MOFid command line and web tools are built on top of other open-source software projects:

* [Open Babel](https://github.com/openbabel/openbabel) cheminformatics toolkit
* [eigen](http://eigen.tuxfamily.org/) is bundled as a dependency for Open Babel
* Make, [cmake](https://cmake.org/overview/), [Node.js](https://nodejs.org/en/), and [Emscripten](https://emscripten.org/index.html) provide the build infrastructure
* [Systre](http://www.gavrog.org/) (and [webGavrog](https://github.com/odf/webGavrog) in the online tool) analyze crystal graph data to assign [RCSR topology symbols](http://rcsr.anu.edu.au/) for MOFs
* [NGL Viewer](https://github.com/arose/ngl) is used to visualize MOF structures and components on the website
* [Kekule.js](http://partridgejiang.github.io/Kekule.js/) enables users to draw molecule substructure queries in the searchdb web tool
