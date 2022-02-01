# Singularity

## Installation
If you have access to a Linux system or high-performance computing cluster, it may be possible to run the MOFid code via [singularity](https://apptainer.org/user-docs/master/quick_start.html), which packages the MOFid installer into a portable, reproducible environment. To get started, refer to documentation from your university or computing center ([example](https://kb.northwestern.edu/page.php?id=85614)) for help on singularity. There may be setup instructions specific to your compute environment. (For example, you may need to load modules or bind paths to set up singularity.)

1. Download the pre-compiled singularity container (`mofid.sif`) from the [most recent release](https://github.com/snurr-group/mofid/releases).
2. Test your installation using `singularity test mofid.sif`. Your installation is successful if you receive the message `Results: 0 errors in 28 MOFs` at the end of the run.

See [additional details](containers.md) about alternate container-based installation methods, such as using [Docker](https://www.docker.com/resources/what-container).

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

MOFID_SIF = "path_to_mofid.sif"
MOF_CIF_TO_ANALYZE = "path_to_mof.cif"
mofid_cmd = ["singularity", "run", MOFID_SIF, "file", MOF_CIF_TO_ANALYZE]
mofid_run = subprocess.run(mofid_cmd, universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
sys.stderr.write(mofid_run.stderr)  # Re-forwarding C++ errors
mofid_output = json.loads(mofid_run.stdout)
```

The `mofid_output` variable above is a dictionary containing eight entries: the MOFid (`mofid`), MOFkey (`mofkey`), SMILES string (`smiles`, `smiles_nodes`, or `smiles_linkers`), topology (`topology`), catenation (`cat`), and basename of the CIF (`cifname`).

## Singularity

See [here](singularity.md) for a quick-start on singularity.

### Installation

If installing singularity on your own Linux-based system, see the [singularity quick start](https://apptainer.org/user-docs/master/quick_start.html); otherwise, contact your system administrator.

The main README uses a prebuilt singularity image (mofid.cif) distributed as a GitHub release. To build it yourself, clone the MOFid repo, then run:

```{bash}
cd path_to_mofid_repo
sudo singularity build ../path_to_output/mofid.sif mofid.def

# singularity build requires root access or remote building
# The build step will also run all the tests

# Once the container is built, you can get a shell into the singularity environment via:
singularity shell mofid.sif
```

### Usage

See the main README for usage instructions.

Alternatively, singularity has a python module called `spython` to spin up containers. Behind the scenes, `spython` uses the `subprocess` package to interface with singularity, but it merges output from stdout and stderr (at the time of writing). The approach in the main README introduces less package dependencies by directly calling `subprocess`. Regardless, the experiments with `spython` are documented below:

```{python}
from spython.main import Client  # installed via: pip install spython
import json

mofid_instance = Client.instance("mofid.sif")
# depending on your HPC environment, you may need to add options=["--bind", "/project_dir:/project_dir"]

# Analyzing a single file
mofid_command = ["python", "/mofid/Python/run_mofid.py", "path_to_mof.cif", "json"]
mofid_stdout_and_stderr = Client.execute(mofid_instance, mofid_command)

# Analyzing a folder is possible, but the output is not as clean

# Clean up
Client.instance_stopall()
```
