#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Load cheminformatics libraries and helper utilities

Loads Open Babel, pybel, and applicable normalization wrappers

@author: Ben Bucior
"""

from easyprocess import EasyProcess  # Install with pip or from https://github.com/ponty/EasyProcess
import sys, os

def path_to_resource(resource):
	# Get the path to resources, such as the MOF DB's or C++ code, without resorting to hardcoded paths
	# However, some absolute paths are still present in extract_moffles.py since they're system-wide
	python_path = os.path.dirname(__file__)
	return os.path.join(python_path, resource)

# Set up local Open Babel data environment before importing the libraries.
# CIF and other SMILES work is handled by the bin/sbu binary, called as a subprocess
os.environ["BABEL_DATADIR"] = path_to_resource("../src/ob_datadir")  # directory with native EOL
TSFM_BIN = path_to_resource("../bin/tsfm_smiles")
OBABEL_BIN = path_to_resource("../openbabel/build/bin/obabel")

def quote(smiles_str):
	# Wraps SMILES strings for Open Babel within single (non-parseable) quotes
	return "'" + smiles_str + "'"

def in_smiles(smiles_str):
	# Adds necessary quotes and Open Babel input notation for SMILES strings
	return "-:" + quote(smiles_str)

def ob_normalize(smiles):
	# Normalizes an arbitrary SMILES string with the same format and parameters as sbu.cpp
	cpp_run = EasyProcess([OBABEL_BIN, in_smiles(smiles), "-xi", "-ocan"]).call()
	cpp_output = cpp_run.stdout
	if (cpp_run.stderr != "1 molecule converted"):
		sys.stderr.write(cpp_run.stderr + "\n")  # Re-fowarding obabel errors
	return cpp_output.rstrip()

def openbabel_replace(mol_smiles, query, replacement):
	# Perform Open Babel transforms, deletions, and/or replacements on a SMILES molecule.
	# With help from on http://baoilleach.blogspot.com/2012/08/transforming-molecules-intowellother.html
	# See also the [Daylight manual on SMARTS](http://www.daylight.com/dayhtml/doc/theory/theory.smarts.html)
	# and phmodel.cpp:208, which clarifies the possibilities of Open Babel replacements
	cpp_run = EasyProcess([TSFM_BIN, quote(mol_smiles), quote(query), quote(replacement)]).call()
	cpp_output = cpp_run.stdout
	sys.stderr.write(cpp_run.stderr)  # Re-fowarding C++ errors
	return ob_normalize(cpp_output.rstrip())

def openbabel_contains(mol_smiles, query):
	# Checks if a molecule (including multi-fragment contains a SMARTS match
	cpp_run = EasyProcess([OBABEL_BIN, in_smiles(mol_smiles), "-s", quote(query), "-xi", "-ocan"]).call()
	cpp_output = cpp_run.stdout
	if (cpp_run.stderr == "1 molecule converted"):
		return True
	elif (cpp_run.stderr == "0 molecules converted"):
		return False
	else:
		sys.stderr.write(cpp_run.stderr + "\n")  # Re-fowarding obabel errors
		return False

