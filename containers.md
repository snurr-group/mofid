# How to run MOFid

This document provides additional information about alternative installation methods for the MOFid package. Some of these commands may require adjustments to run on your system.


## Docker

The initial containerization experiments for MOFid used [Docker](https://www.docker.com/resources/what-container). Installation instructions for Docker are highly system-dependent and beyond the scope of this document. For example, Windows users may need to set up WSL2 or Virtualbox, depending on software licenses and versions. Two potentially helpful links for Linux include the [installation](https://docs.docker.com/engine/install/ubuntu/) and [configuration](https://docs.docker.com/engine/security/rootless/) docs.

### Installation

1. Once Docker is set up, clone the MOFid repository.
2. Build the container by running `docker build -t mofid path_to_mofid_repo`

### Usage

To start a shell within the Docker container, run `docker run -it mofid /bin/bash`

To analyze a single MOF crystal structure:

```{bash}
# Configure paths in the first two lines, then run this block
MOFID_IN_FILENAME="P1-Cu-BTC.cif"
MOFID_IN_DIR="$PWD"  # absolute path to parent directory, e.g. current directory via $PWD

docker run \
    --mount type=bind,source="$MOFID_IN_DIR",target="$MOFID_IN_DIR",readonly \
    mofid \
    python /mofid/Python/run_mofid.py "$MOFID_IN_DIR/$MOFID_IN_FILENAME" /mofid/TempOutput json
```

Or, for analyzing an entire folder:

```{bash}
# Configure the paths in the first two lines, then run this block
MOFID_IN_DIR="$PWD/Resources/TestCIFs"
MOFID_OUT_DIR="$PWD/mofid_output"

mkdir -p "$MOFID_OUT_DIR"
docker run \
    --mount type=bind,source="$MOFID_IN_DIR",target=/data,readonly \
    --mount type=bind,source="$MOFID_OUT_DIR",target=/out \
    mofid \
    Scripts/run_folder.sh /data /out
```


## Singularity

See the [main project README](https://github.com/snurr-group/mofid/blob/master/containers.md) for a quick start on singularity.

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


## Install on base system

This software can also be installed directly as a Python package without using Docker, singularity, or other containerization packages. To install and run mofid, see the instructions and installation requirements below. (from the README in mofid 1.0.1)

### Requirements
1. A Python environment is required. If you do not have a Python environment installed, we recommend downloading and installing [Anaconda](https://www.anaconda.com/distribution/#download-section). MOFid is compatible with both Python 2/3.
2. Make sure you have the following: a C++ compiler, [`cmake`](https://cmake.org/), and access to GNU commands (such as `make`). These are all typically available on Linux machines. If running on Windows, we recommend using [Cygwin](https://www.cygwin.com/) and including the `cmake`, `make`, `wget`, `gcc-core`, `gcc-g++`, and `pkg-config` packages in addition to the default options during the Cygwin installation process.
3. Make sure you have the [Java Runtime Environment](https://www.java.com/en/download/) installed and included in your system's path. If unsure, try running `java` in the command line to see if it successfully calls Java.

### Installation
1. Run `make init` in the base `/mofid` directory.
2. Run `python set_paths.py` followed by `pip install .` in the base `/mofid` directory.  If you encounter permissions errors (typically not with Anaconda), you may need to run `pip install --user .`

### Usage
In a Python script, the user simply has to call the `run_mofid.cif2mofid(cif_path,output_path='Output')` function. The first argument is required and is the path to the MOF CIF. The second argument is optional and is the directory to store the MOFid decomposition information, which defaults to `/Output` if not specified. An example of how to call MOFid is shown below.

```python
from mofid.run_mofid import cif2mofid
cif_path = '/path/to/my/mof.cif'
mofid = cif2mofid(cif_path)
```

The output of the `mofid.cif2mofid` function is a dictionary containing eight entries: the MOFid (`mofid`), MOFkey (`mofkey`), SMILES string (`smiles`, `smiles_nodes`, or `smiles_linkers`), topology (`topology`), catenation (`cat`), and basename of the CIF (`cifname`).

