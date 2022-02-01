# MOFid
A system for rapid identification and analysis of metal-organic frameworks.

Please cite [DOI: 10.1021/acs.cgd.9b01050](https://pubs.acs.org/doi/abs/10.1021/acs.cgd.9b01050) if you use MOFid in your work.

## Objective
Supplement the current MOF naming conventions with a canonical, machine-readable identifier to facilitate data mining and searches. Accomplish this goal by representing MOFs according to their nodes + linkers + topology

## Usage and Installation Instructions
There are three main ways in which you can use MOFid:
1. From your browser.
2. By compiling the MOFid source code and running it locally.
3. By using Singularity or Docker to run a pre-built image of the MOFid code locally.

### Browser-Based MOFid
Visit [https://snurr-group.github.io/web-mofid](https://snurr-group.github.io/web-mofid/) to quickly and easily run MOFid in your browser! No programming skills are required.

### Compiling from Source
See [compiling.md](compiling.md) for how to compile and run MOFid from source.

### Containerized MOFid
See [singularity.md](singularity.md) for how to run MOFid via a Singularity container.

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
