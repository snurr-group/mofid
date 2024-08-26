---
layout: default
title: Singularity
permalink: /singularity/
---

# Singularity
{: .no_toc}

## Table of Contents
{: .no_toc .text-delta}

1. TOC
{:toc}

## Installation
If you have access to a Linux machine or HPC cluster, it may be possible to run MOFid via [Singularity](https://apptainer.org/user-docs/master/quick_start.html), which packages the MOFid installer into a portable and reproducible environment. To get started, refer to documentation from your university or computing center ([example](https://kb.northwestern.edu/page.php?id=85614)) for help on singularity. There may be setup instructions specific to your compute environment. E.g. you may need to load modules or bind paths to set up Singularity.

1. Download the pre-compiled Singularity container `mofid.sif` from the [most recent release](https://github.com/snurr-group/mofid/releases).
2. Test your installation using `singularity test mofid.sif`. Your installation is successful if you receive the message `Results: 0 errors in 28 MOFs`.

Alternatively, if you wish to use Docker see [docker.md](docker.md).

## Usage
The Singularity container wraps all of MOFid into a single package.

As a command line tool...

```bash
# Analyzing a single MOF crystal structure
./mofid.sif file path_to_mof.cif
# alternatively: singularity run mofid.sif file path_to_mof.cif

# Analyzing a folder
./mofid.sif folder path_to_input_cif_folder path_to_mofid_output
# By default, path_to_mofid_output is set to "Output/" in your current directory
```

or as part of a Python script...

```python
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
