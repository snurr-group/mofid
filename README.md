# MOFid
A system for rapid identification and analysis of metal-organic frameworks

## Objective
**Supplement** the current MOF naming conventions with a canonical, machine-readable identifier to facilitate data mining and searches.  Accomplish this goal by representing MOFs according to their nodes + linkers + topology

## Installation notes/requirements
* Needs a C++ compiler, make, cmake, and a Java Runtime Environment (for Systre)
* sqlite3 is useful for analysis.
* Testing code is based on Python 2.7.  Long-term it would be good to move everything to Python 3.x, but for now Open Babel is primarily a Python 2 library so everything else is currently kept that way.  Based on the features used by the code, the biggest compatibility change in the future will likely be `print` statements.
* The Python test suite requires the `easyprocess` package to call executables like MOFid or Systre.  Depending on your environment, this can be installed using a command like `pip install --user easyprocess`
* Systre is required for topology detection (check that the path is correct in Python/extract_moffles.py).  `jq` is a useful utility to parse json files but is not strictly necessary.  Both of these can be automatically set up using `make download`
* The Python code directly calls a C++ executable in most cases instead of directly using the `openbabel` Python library.  If `openbabel` or `pybel` is needed but raise an error, you may need to modify the library path, e.g. `export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH`
* To get started, run `make init` the first time to compile Open Babel and set up the cmake environment.  Modifications to `sbu.cpp` or other code can be quickly compiled using `make bin/sbu` and tested using `make test`.

## TODO
* Add an intro to the project structure and tutorial on basic usage.
* Include a link to the paper for additional documentation.


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



