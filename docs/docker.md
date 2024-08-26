---
layout: default
title: Docker
permalink: /docker/
---

# Docker
{: .no_toc}

## Table of Contents
{: .no_toc .text-delta}

1. TOC
{:toc}

## Installation

Installation instructions for [Docker](https://www.docker.com/) are beyond the scope of this document. Two potentially helpful links for Linux include the [installation](https://docs.docker.com/engine/install/ubuntu/) and [configuration](https://docs.docker.com/engine/security/rootless/) pages.


1. Once Docker is set up, clone the MOFid repository.
2. Build the container by running `docker build -t mofid path_to_mofid_repo`

## Usage

To start a shell within the Docker container, run `docker run -it mofid /bin/bash`

To analyze a single MOF crystal structure...

```bash
# Configure paths in the first two lines, then run this block
MOFID_IN_FILENAME="P1-Cu-BTC.cif"
MOFID_IN_DIR="$PWD" # absolute path to parent directory, e.g. current directory via $PWD

docker run \
    --mount type=bind,source="$MOFID_IN_DIR",target="$MOFID_IN_DIR",readonly \
    mofid \
    python /mofid/Python/run_mofid.py "$MOFID_IN_DIR/$MOFID_IN_FILENAME" /mofid/TempOutput json
```

or to analyze an entire folder...

```bash
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
