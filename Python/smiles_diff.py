#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Run a diff between two SMILES codes

Report common classes of errors in the calculated MOFid.

@author: Ben Bucior
"""

import sys, os

# TODO: Refactor the OpenBabel loading as another helper import
def path_to_resource(resource):
	# Get the path to resources, such as the MOF DB's or C++ code, without resorting to hardcoded paths
	# However, some absolute paths are still present in extract_moffles.py since they're system-wide
	python_path = os.path.dirname(__file__)
	return os.path.join(python_path, resource)
os.environ["BABEL_DATADIR"] = path_to_resource("../src/ob_datadir")
import pybel  # Read SMILES to calculate molecular formulas, etc.


def multi_smiles_diff(smiles1, smiles2):
	# What are the differences between two dot-separated SMILES codes, a reference and candidate?
	# Just bond orders?  The entire bond structure?  Nothing similar?

	if smiles1 == smiles2:
		return []
	else:
		return ["smiles"]  # STUB: Let's focus on single diff's, then come back to the more complicated case

	# TODO: complete this function
	diffs = []
	# Consider making the results from single_smiles_diff as an enum, class, or other rankable level of error for
	# prioritizing which "match" is closest

	parts1 = smiles1.split('.')
	parts2 = smiles2.split('.')

	categorized = []
	for code1 in parts1:
		for code2 in parts2:
			if code1 == code2:
				categorized.append(["equal", code1, code2])
	for match in categorized:
		parts1.remove(match[1])
		parts2.remove(match[2])

	print categorized
	print parts1, parts2

	processed1 = [x[1] for x in categorized]
	processed2 = [x[2] for x in categorized]



	# LOOP OVER THINGS LEFT OVER
	# First diff the two codes to see if it eliminates one of smiles2
	# If not, extract the opposite node from "categorized" and see if it belongs to that (fallback?)
	# If that fails, report an "extra" smiles component


	diffs.append("smiles")  # Maybe something about extracting all the categorized[x][0] terms?
	return diffs


def single_smiles_diff(smiles1, smiles2):
	# What are the differences between two single-component SMILES?
	if '.' in smiles1 or '.' in smiles2:
		raise ValueError("Only a single component is allowed")
	if smiles1 == smiles2:
		return "none"
	error_codes = ["ERROR", "NA", ""]
	if smiles1 in error_codes or smiles2 in error_codes:
		return "ERROR"

	mol1 = pybel.readstring("smi", smiles1)
	mol2 = pybel.readstring("smi", smiles2)

	if mol1.formula != mol2.formula:
		return "formula"

	def is_organic(mol):
		return 'C' in mol.OBMol.GetSpacedFormula().split(' ')
	if is_organic(mol1) and is_organic(mol2):
		molec_type = "linker"
	else:
		molec_type = "node"

	# Molecules with the same single bond structure appear to have the same canonical atom order, just
	# different bond orders, bracketing, and explicit hydrogen atoms.
	base1 = smiles1
	base2 = smiles2
	for special in 'H$[]()@-=#:$/\\':
		base1 = base1.replace(special, '')
		base2 = base2.replace(special, '')
	base1 = base1.upper()
	base2 = base2.upper()

	if base1 != base2:
		return molec_type + "_single_bonds"
	else:
		return molec_type + "_bond_orders"


def usage():
	raise SyntaxError("Run a diff between two SMILES strings.  Be sure to quote to escape bash globbing, etc.")


if __name__ == "__main__":
	args = sys.argv[1:]
	if len(args) != 2:
		usage()

	print single_smiles_diff(args[0], args[1])  # FIXME: use multi_ once implemented
