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


DIFF_LEVELS = dict({
	'same' : 0,
	'formula' : 50,
	'linker_bond_orders' : 10,
	'node_bond_orders' : 11,
	'linker_single_bonds' : 20,
	'node_single_bonds' : 21,
	'ERROR' : 99,
	'ERR_MAX' : 99999
})


def multi_smiles_diff(smiles1, smiles2):
	# What are the differences between two dot-separated SMILES codes, a reference and candidate?
	# Just bond orders?  The entire bond structure?  Nothing similar?
	if smiles1 == smiles2:
		return []

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


	# TODO: refactor inner parts of loop possibly as:
	# find_closest_match(smiles, unmatched_list, extra_list)

	while len(parts2):
		current = parts2.pop()  # MAYBE ENCAPSULATE THE REST AS A SEMI-LOOP, THEN RE-INCORPORATE??
		best_match = ['ERR_MAX', '', True]  # Error level, SMILES, unpaired

		processed1 = [x[1] for x in categorized]
		processed2 = [x[2] for x in categorized]

		for unpaired in parts1:
			test_err = single_smiles_diff(current, unpaired)
			if DIFF_LEVELS[test_err] < DIFF_LEVELS[best_match[0]]:
				best_match = [test_err, unpaired, True]
		for repeat in processed1:
			test_err = single_smiles_diff(current, repeat)
			if DIFF_LEVELS[test_err] < DIFF_LEVELS[best_match[0]]:
				best_match = [test_err, repeat, False]

		if not best_match[2]:
			best_match[0] = best_match[0] + "_extra"  # Duplicated node/linker, likely from inconsistent representations
		if best_match[0] != 'ERR_MAX':
			categorized.append([best_match[0], best_match[1], current])

	# I'm doing parts2 first, but the same process applies to parts1
	# Rewrite as a function, which includes an index of which element to pursue????

	#For a stub, just return ["smiles"]
	return [x[0] for x in categorized if x[0] != 'equal']


def single_smiles_diff(smiles1, smiles2):
	# What are the differences between two single-component SMILES?
	if '.' in smiles1 or '.' in smiles2:
		raise ValueError("Only a single component is allowed")
	if smiles1 == smiles2:
		return "equal"
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

	print multi_smiles_diff(args[0], args[1])
