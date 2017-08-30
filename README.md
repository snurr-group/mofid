# MOFid
Toward a universal MOF identifier

## Objective
**Supplement** the current MOF naming conventions with a canonical, machine-readable identifier to facilitate data mining and searches.  Accomplish this goal by representing MOFs according to their nodes + linkers + topology

## TODO
* Add an intro to the project structure and tutorial on basic usage.
* Give background on validating MOFid code against ToBaCCo and GA MOFs.
* Include paper outline somewhere?
* Clean up repo and comb through old known issues and potentially useful papers for additional validation.

## Known issues (incomplete list)
### Technical questions
* What should the final MOFid format look like?  SMILES has many advantages (interoperability, interpretability) but also downsides (search engine friendliness, not compact, sensitivity to aromaticity models).
* How are counter anions handled?  How should they be classified in the ID?
* Further consideration of how nodes are defined.  Metal oxides?  This question may change the assigned topology and simplification of some structures.

### Missing features
* Only P1 CIFs are recognized.  Symmetry operations are not applied
* Structures are read as-is.  Hydrogens are not (yet) automatically added
* The HTML (emscripten) interface is still a bare-bones proof-of-concept.  Visualization and potential integration with MOF Explorer is still necessary.

### Major bugs
* Handling of formal charges, which is recommended for SMILES and even InChI.  Currently only carboxylates and certain nitrogen rings are handled by SMARTS pattern recognition.
* Additional topological analysis is required.  Existing RCSR topologies are only obtained for ~2,000 of the CIFs from CoRE MOF 2.0, which is a similar fraction as an automated TOPOS analysis.
* A few structures crash sbu.exe or the Python code.
* Look for warnings in the CoRE MOF output.  These aren't inherently bad in a semi-messy data source, but they also flag potentially bad assumptions in the MOFid code.

### Bugs with upstream Open Babel project
* MOFid code uses the Open Babel project, as downloaded from Github on 6/23/16, commit 03fa0761dbd3b10a6d31036753faf8dddeda7ac5.  See also their [C++ tutorial](http://openbabel.org/wiki/Developer:Cpp_Tutorial)
* Steps are in progress to merge changes, such as periodic boundary conditions, back into the upstream project.  Then, this dependency can be imported as a git submodule.  In the meantime, there are a few incompatibilities between MOFid and the upstream Open Babel project.
* Open Babel's aromaticity model has problems with certain nitrogen-containing rings.  Certain heteroaromatic linkers from ToBaCCo are classified with incorrect bond orders.  Resolving this [issue](https://github.com/openbabel/openbabel/issues/1360) is in progress.
* By default, Open Babel's bond assignment code (`OBMol::ConnectTheDots()`) deletes bonds in excess of a "maximum valence."  That means bonds can unexpectedly disappear from some CIFs.  Due to the way a few details are implemented, that means (1) assigning neighbors and topology for some metals is inconsistent, and (2) hydroxyls are not reported correctly due to some other workarounds in MOFid.


## Misc. notes
### Eclipse
Generated the eclipse code hinting by running `cmake -G "Eclipse CDT4 - Unix Makefiles" ../src` from within bin/ based on a [helpful post](http://stackoverflow.com/questions/11645575/importing-a-cmake-project-into-eclipse-cdt).  Another [website](http://www.badprog.com/c-eclipse-installation-of-c-c-development-tools-cdt-and-cygwin-for-windows) might be useful if I try to get the "Run" button working later.

The [main API reference](http://openbabel.org/dev-api/namespaceOpenBabel.shtml) and [OBMol class reference](http://openbabel.org/dev-api/classOpenBabel_1_1OBMol.shtml) are also great resources.  Also see the [installation instructions](https://openbabel.org/docs/dev/Installation/install.html#local-build).  With the new makefile configurations, run `make openbabel/fast` to make a quick local build, and the new sbu Makefile will automatically compile against the local directory.

### OpenBabel listserv
Listserv entries discussing periodic boundary conditions.

* [Periodic boundaries and atom typing large systems (UFF)](https://www.mail-archive.com/openbabel-discuss@lists.sourceforge.net/msg02002.html)
* [OpenBabel and Periodic Systems (multiple molecules per OBMol)](https://sourceforge.net/p/openbabel/mailman/message/7048390/)
* [Crystallography help (fractional-cartesian conversion matrix)](https://sourceforge.net/p/openbabel/mailman/message/7049196/)

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



