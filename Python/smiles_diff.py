#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Run a diff between two SMILES codes

Report common classes of errors in the calculated MOFid.

@author: Ben Bucior
"""

import sys, os
import re

from cpp_cheminformatics import ob_normalize, openbabel_replace, openbabel_formula, openbabel_GetSpacedFormula


DIFF_LEVELS = dict({
	'same' : 0,
	'formula' : 50,
	'nonplanar_carboxylate': 5,
	'phenyl_radicals': 6,
	'fg_bond_location' : 8,
	'linker_bond_orders' : 10,
	'node_bond_orders' : 11,
	'linker_single_bonds' : 20,
	'node_single_bonds' : 21,
	'charges' : 40,
	'order' : 90,
	'ERROR' : 99,
	'ERR_MAX' : 99999
})


def find_closest_match(smiles, preferred_list, extra_list):
	# Find the closest match for a SMILES in two candidate lists.
	# For each candidate, run single_smiles_diff, and rank the
	# possible errors by their value in DIFF_LEVELS.
	best_match = ['ERR_MAX', '', True]  # [Error level, SMILES, is this from the extra list?]

	for option in preferred_list:
		test_err = single_smiles_diff(smiles, option)
		if DIFF_LEVELS[test_err] < DIFF_LEVELS[best_match[0]]:
			best_match = [test_err, option, False]
	for repeat in extra_list:
		test_err = single_smiles_diff(smiles, repeat)
		if DIFF_LEVELS[test_err] < DIFF_LEVELS[best_match[0]]:
			best_match = [test_err, repeat, True]

	if best_match[2]:
		best_match[0] = best_match[0] + "_extra"  # Duplicated node/linker, likely from inconsistent representations
	if best_match[0] == 'ERR_MAX_extra':
		best_match = None

	return best_match


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

	while len(parts1):
		current = parts1.pop()
		processed2 = [x[2] for x in categorized]
		best_match = find_closest_match(current, parts2, processed2)
		if best_match is not None:
			categorized.append([best_match[0], current, best_match[1]])
			if not best_match[2]:  # If part of the "preferred list" (parts2)
				parts2.remove(best_match[1])

	while len(parts2):
		current = parts2.pop()
		processed1 = [x[1] for x in categorized]
		best_match = find_closest_match(current, parts1, processed1)
		if best_match is not None:
			categorized.append([best_match[0], best_match[1], current])
			assert best_match[2]  # parts1 is empty, so the best match must have already been added to categorized

	#For a stub, just return ["smiles"]
	err_codes = []
	for x in categorized:
		err = x[0]
		if err != 'equal' and err not in err_codes:
			err_codes.append(err)

	if len(err_codes) == 0 and categorized[0][0] == 'equal' and smiles1 != smiles2:
		return ['order']

	return err_codes


def single_smiles_diff(smiles1, smiles2):
	# What are the differences between two single-component SMILES?
	if '.' in smiles1 or '.' in smiles2:
		raise ValueError("Only a single component is allowed")
	if smiles1 == smiles2:
		return "equal"
	error_codes = ["ERROR", "NA", "", "*"]
	if smiles1 in error_codes or smiles2 in error_codes:
		return "ERROR"

	mol1_formula = openbabel_formula(smiles1)
	mol2_formula = openbabel_formula(smiles2)

	if re.sub('[+-]', '', smiles1) == re.sub('[+-]', '', smiles2):
		return "charges"

	def strip_extra(formula):
		# Strip hydrogens and charges from molecular formula.
		# We don't have to worry about greatest common factor, etc., since it's absolute atom counts.
		return re.sub('[+-]', '', re.sub(r'H\d+', '', formula))
	if strip_extra(mol1_formula) != strip_extra(mol2_formula):
		return "formula"

	def radical_to_carb(smiles):
		return openbabel_replace(smiles, '[#6:1][#6D3:2](~[O-0:3])[O:4]', '[#6:1][#6:2](=[O:3])[O-:4]')
	if radical_to_carb(smiles1) == radical_to_carb(smiles2):
		return "nonplanar_carboxylate"

	if ob_normalize(smiles1.replace('[c]', 'c')) == ob_normalize(smiles2.replace('[c]', 'c')):
		return "phenyl_radicals"

	def move_hydrogen(smiles):
		# Transfer a proton from a carbonyl to a nearby aromatic carbon (and/or another linker)
		# Shows up in certain linkers when functional groups are assigned to the wrong neighbor
		return ob_normalize(smiles.replace('[OH]', 'O').replace('[c]', 'c'))
	if move_hydrogen(smiles1) == move_hydrogen(smiles2):
		return "fg_bond_location"

	def is_organic(smiles):
		return 'C' in openbabel_GetSpacedFormula(smiles).split(' ')
	if is_organic(smiles1) and is_organic(smiles2):
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
