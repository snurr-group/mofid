# MOFid
A system for rapid identification and analysis of metal-organic frameworks

**WARNING: code is under active development and refactoring.  Details subject to change**

## Objective
**Supplement** the current MOF naming conventions with a canonical, machine-readable identifier to facilitate data mining and searches.  Accomplish this goal by representing MOFs according to their nodes + linkers + topology

## Installation notes/requirements
1. Make sure you have the following: a C++ compiler installed, [cmake](https://cmake.org/), and access to GNU commands (such as `make`) . If running a Windows machine, this can all be done by installing [Cygwin](https://www.cygwin.com/), making sure to include both the `make` and `wget` packages in addition to the default options.
2. Make sure you have a [Java Runtime Environment](https://www.java.com/en/download/) installed, which is needed for [Systre](http://gavrog.org/).
Java in path

3. Install [Systre](https://github.com/odf/gavrog/releases), which is required for topology detection (check that the path is correct in the `/Python/extract_moffles.py` file). Both of these can be automatically downloaded and installed using `make download`.
4. Run `make init` to compile a customized development version of Open Babel and set up the `cmake` environment.  Modifications to `sbu.cpp` or other codes can be quickly compiled using `make exe` and tested using `make test`. `TODO: Why would one need this (the last sentence)?`

Additional notes:
1. The Python code directly calls a C++ executable in most cases instead of directly using the `openbabel` Python library.  If `openbabel` or `pybel` is needed but raise an error, you may need to modify the library path, e.g. `export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH`
2. [SQlite3](https://www.sqlite.org/index.html) can be useful useful for analysis.
3. The Python test suite requires the [`subprocess32`](https://pypi.org/project/subprocess32/) package to call executables like MOFid or Systre, which also backports functionality like call timeouts from Python 3.x.  Depending on your environment, this can be installed using a command like `pip install subprocess32==3.5.3`. `TODO: Make Py3 compatible and remove dependence on subprocess32.`

## Code overview
* openbabel/ contains a modified version of the upstream openbabel library (development version)
* src/ is the MOFid source code.
	* sbu.cpp is the main MOFid deconstruction.  By default, it writes relevant CIFs and a Systre topology.cgd file to an Output/ directory
	* sobgrep.cpp does SMARTS substructure searches on a single bonded version of the MOF
	* searchdb.cpp is a utility to search through a SMILES database file
	* tsfm_smiles.cpp applies Open Babel's version of reaction transformations to a target SMILES molecule.  It was developed primarily to remove the openbabel.py dependency for the Python source code and avoid incompatibility issues between different Open Babel versions (particularly for atoms with nonstandard valence, like the oxygen in the MOF-5 node)
	* Other files are helpful classes so that sbu.cpp is no longer ~2k lines of code
* Python scripts are commonly used for validation and chaining together sbu and systre outputs
* Bash scripts are misc. utilities, including analyzing databases on high-throughput supercomputers
* The R code is currently unused but may be revisited when writing the manuscript.


## TODO
* Include a link to the paper for additional documentation.
* Later: will include notes on how to compile the web version and host it on Github.  Will be revisited after finishing more urgent analysis and code.


## Misc. notes
### Open Babel API notes
See the [C++ tutorial](http://openbabel.org/wiki/Developer:Cpp_Tutorial) for Open Babel.  Their [main API reference](http://openbabel.org/dev-api/namespaceOpenBabel.shtml) and [OBMol class reference](http://openbabel.org/dev-api/classOpenBabel_1_1OBMol.shtml) are also great resources.  Also see the [installation instructions](https://openbabel.org/docs/dev/Installation/install.html#local-build).  With the new makefile configurations, run `make openbabel/fast` to make a quick local build, and the new sbu Makefile will automatically compile against the local directory.

### Emscripten installation notes
In Windows, download the [portable Emscripten SDK](http://kripken.github.io/emscripten-site/docs/getting_started/downloads.html#platform-notes-installation-instructions-portable-sdk) and install using Git Bash.  The general process is described on the website:

```
Follow instructions:
# Fetch the latest registry of available tools.
./emsdk update

# Download and install the latest SDK tools.
./emsdk install latest

# Make the "latest" SDK "active"
./emsdk activate latest
```

Setting the paths doesn't quite work due to mismatches between CMD and bash, but a new [importer script](Scripts/import_emscripten.sh) will automatically take care of those details in Git Bash, if sourced.  Sourcing is critical to ensure that the changes to environment variables propagate to the current shell instead of being confined to a subprocess.  Then, the Emscripten `make` processes are automated inside the `web` and `init-web` targets of the project Makefile.

Additional emscripten links and diagnostics (like `emcc -v`) are sketched in the corresponding [notebook directory](Notebooks/20170810-emscripten/emscripten_installation.txt).



