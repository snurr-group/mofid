# MOFid
A system for rapid identification and analysis of metal-organic frameworks.

Check out [https://snurr-group.github.io/web-mofid/](https://snurr-group.github.io/web-mofid/) to quickly and easily run these tools in your browser! No programming skills required. If you wish to generate the MOFid for a large number of structures, see the Python-based interface described below.

Please cite [DOI: 10.1021/acs.cgd.9b01050](https://pubs.acs.org/doi/abs/10.1021/acs.cgd.9b01050) if you use MOFid in your work.

## Objective
Supplement the current MOF naming conventions with a canonical, machine-readable identifier to facilitate data mining and searches. Accomplish this goal by representing MOFs according to their nodes + linkers + topology

## Installation
If you have access to a Linux system or high performance computing cluster, it may be possible to run the MOFid code via [singularity](https://apptainer.org/user-docs/master/quick_start.html), which packages the mofid installer into a portable, reproducible environment. To get started, refer to documentation from your university or computing center ([example](https://kb.northwestern.edu/page.php?id=85614)) for help on singularity. There may be setup instructions specific to your compute environment. (For example, you may need to load modules or bind paths to set up singularity.)

1. Set up singularity and test your installation, e.g. `singularity exec library://ubuntu cat /etc/lsb-release`
2. Download the pre-compiled singularity container from GitHub, e.g. via `singularity pull mofid.sif GITHUB_URL_TODO`
3. Test your installation using `singularity test mofid.sif`

See [additional details](https://github.com/snurr-group/mofid/blob/master/containers.md) about alternate installation methods, such as using [Docker](https://www.docker.com/resources/what-container) or compiling the Python package yourself.

## Usage
The singularity container wraps all of the MOFid software into a single package.

As a command line tool:

```{bash}
# Analyzing a single MOF crystal structure
./mofid.sif file path_to_mof.cif
# alternatively: singularity run mofid.sif file path_to_mof.cif

# Analyzing a folder
./mofid.sif folder path_to_input_cif_folder path_to_mofid_output
# By default, path_to_mofid_output is set to "Output/" in your current directory
```

Or, as part of a Python script:

```{python}
import json
import sys

import subprocess
# If using Python 2, you may need to install subprocess32 and import via:
# import subprocess32 as subprocess

mofid_cmd = ["singularity", "run", "python", "/mofid/Python/mofid_json.py", "path_to_mof.cif"]
mofid_run = subprocess.run(mofid_cmd, universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
sys.stderr.write(mofid_run.stderr)  # Re-forwarding C++ errors
mofid_output = json.loads(mofid_run.stdout)
```

The `mofid_output` variable above is a dictionary containing eight entries: the MOFid (`mofid`), MOFkey (`mofkey`), SMILES string (`smiles`, `smiles_nodes`, or `smiles_linkers`), topology (`topology`), catenation (`cat`), and basename of the CIF (`cifname`).

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
