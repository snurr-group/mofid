Notebook for MOF barcode system
========

Objective
---------
Decompose a MOF into its nodes, linkers, and topology to perform additional analysis (comparisons, chemoinformatics on its components, easy topology, etc.)


TODO
----
* OpenBabel does not currently handle PBC in bond perception.  Would need to make a new ConnectTheBonds with periodic boundaries for CIFs (see mol.cpp:3534).  Also to keep the oxygen with tetrahedral bonding, etc, probably as an option to the CIF parser.


Ideas
-----
* See also <http://drivendata.github.io/cookiecutter-data-science> for project organization ideas


Notes
-----
Generated the eclipse code hinting by running `cmake -G "Eclipse CDT4 - Unix Makefiles" ../src` from within bin/ based on a [helpful post](http://stackoverflow.com/questions/11645575/importing-a-cmake-project-into-eclipse-cdt).  Another [website](http://www.badprog.com/c-eclipse-installation-of-c-c-development-tools-cdt-and-cygwin-for-windows) might be useful if I try to get the "Run" button working later.




Table of contents
-----------------
* openbabel: Main openbabel distribution from Github on 6/23/16, commit 03fa0761dbd3b10a6d31036753faf8dddeda7ac5.  See <http://openbabel.org/wiki/Install_(Cygwin)> for installation instructions and <http://openbabel.org/wiki/Developer:Cpp_Tutorial> for the C++ tutorial.


Reading
-------
* 

