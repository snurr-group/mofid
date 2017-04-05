#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Calculate MOF linkers

Use the CSD criteria (no bonds to metals) in my modified OpenBabel code and
sbu.cpp to decompose MOFs into fragments.  Compare actual fragments from
ToBACCo or hMOF structures against their "recipe."

@author: Ben Bucior
"""

import os, sys
import subprocess
import glob
import json
# import copy  # copy.deepcopy(x)
# import re
import openbabel  # for visualization only, since my changes aren't backported to the python library
from extract_moffles import cif2moffles, assemble_moffles, parse_moffles, compare_moffles, extract_linkers

SBU_BIN = "C:/Users/Benjamin/Git/mofid/bin/sbu.exe"
HMOF_DB = "C:/Users/Benjamin/Git/mofid/Resources/hmof_linker_info.json"
TOBACCO_DB = "C:/Users/Benjamin/Git/mofid/Resources/tobacco_info.json"
KNOWN_DB = "C:/Users/Benjamin/Git/mofid/Resources/known_info.json"
DEFAULT_CIFS = "Data/tobacco_L_12/quick"


def any(member, list):
	# Is the member any part of the list?
	return member in list

def summarize(results):
	# Summarize the error classes for MOFFLES results
	summarized = {'mofs': results, 'errors': dict()}
	error_types = {'topology': 0, 'smiles': 0, 'both': 0, 'success': 0, 'undefined': 0}
	for match in results:
		if match['match'] == 'NA':  # NA cases where we don't know the nodes and linkers
			assert 'Undefined composition' in match['errors']
			error_types['undefined'] += 1
		elif match['match']:
			error_types['success'] += 1
		elif len(match['errors']) > 1:
			error_types['both'] += 1
		elif 'topology' in match['errors']:
			error_types['topology'] += 1
		elif 'linkers' in match['errors']:
			error_types['smiles'] += 1
	summarized['errors']['error_types'] = error_types
	summarized['errors']['total_cifs'] = len(results)
	return summarized

class MOFCompare:
	# Compares a MOF's decomposition against its expected components
	def __init__(self):
		self.db_file = None
		# self.load_components()  # Used in subclasses
		raise UserWarning("MOF parsing only (partially) implemented for hMOFs and ToBACCo")

	def load_components(self, db_file = None):
		# Load node and linker compositions from a saved JSON definition (to keep this file cleaner)
		if db_file is None:
			db_file = self.db_file
		with open(db_file, "r") as inp:
			mof_db = json.load(inp)
		self.mof_db = mof_db
		return None

	def extract_db_smiles(self, db_dict):
		# Save all the nodes and linker definitions to the current directory
		# Then, convert to a big svg with
		# obabel linkers.can -O linkers.svg -xe -xC
		# Note: this function and the json file will be deprecated by the metal+linker
		# split paradigm suggested at group meeting ("molecule subtraction")
		linkers = db_dict["linkers"]
		for id in linkers:
			print(linkers[id] + " " + id)


class KnownMOFs(MOFCompare):
	def __init__(self):
		pass

	def parse_filename(self, mof_path):
		# Extract basename of the MOF.  test_cif will convert it to the reference moffles string
		return os.path.splitext(os.path.basename(mof_path))[0]  # Get the basename without file extension

	# def expected_moffles?  test_cif is used in ToBACCo.


class HypoMOFs(MOFCompare):
	def __init__(self):
		self.db_file = HMOF_DB
		raise ValueError("hMOF comparison not yet implemented.  Stub code needs to be updated.")
		self.mof_db = None
		self.load_components()

	def parse_filename(self, hmof_path):
		# Extract hMOF recipes from the filename, formatted as xxxhypotheticalMOF_####_i_#_j_#_k_#_m_#.cif
		codes = {"num": None, "i": None, "j": None, "k": None, "m": None, "cat": 0}
		parts = hmof_path.strip(".cifCIF").split("_")
		flag = None
		for part in parts:
			if flag is not None:
				codes[flag] = part  # For now, keep it as a string, since our json will be string-based, anyway
				flag = None
				continue
			if part.count("MOF") > 0 or part.count("mof") > 0:
				flag = "num"
				continue
			if part in codes.keys():
				flag = part
				continue
		return codes

	def stub(self):
		# Old code about hMOFs
		###hmof_db = load_components(SBU_DB)
		# print hmof_db["nodes"]["0"][1]
		new_linkers = []
		CIF_DIR = "C:/Users/Benjamin/Desktop/Jiayi/Files/Dataset Comparison/hMOF"
		for cif_file in glob.glob(CIF_DIR + '/*.[Cc][Ii][Ff]'):
			# Handle with hMOF class: id = parse_filename(cif_file)
			### STUB:
			id = dict()
			if id['m'] != "0" or id["i"] != "0" or id["k"] != id["j"] or id["j"] in new_linkers:
				# Later, this will probably just be id['m']!='0' since we'll want to temporarily exclude functionalization
				continue
			print cif_file
			for fragment in extract_linkers(cif_file):
				print fragment
			new_linkers.append(id["j"])
			print "\n"


class TobaccoMOFs(MOFCompare):
	def __init__(self):
		self.db_file = TOBACCO_DB
		self.load_components()

	def parse_filename(self, tobacco_path):
		# Extract ToBACCo recipes from the filename
		# Format: topology_sym_x_node_type_sym_x_node2_type_L_linkernum.cif ("_" for empty linker)
		# Can we parse this using ToBACCo's own code??
		codes = {"name": None, "nodes": [], "linker": None, "topology": None}
		mof_info = os.path.splitext(os.path.basename(tobacco_path))[0]  # Get the basename without file extension
		codes['name'] = mof_info

		parsed = mof_info.split("_", 1)
		codes['topology'] = parsed[0]
		mof_info = parsed[1]

		while "sym_" in mof_info:
			parsed = mof_info.split("_")
			codes['nodes'].append("_".join(parsed[0:4]))
			mof_info = "_".join(parsed[4:])

		# Not sure why bcs, etc., have an extra underscore in the topology.
		mof_info = mof_info.strip("_")
		if mof_info == "":
			mof_info = "L__"
		codes['linker'] = mof_info

		return codes

	def expected_moffles(self, cif_path):
		# What is the expected MOFFLES based on the information in a MOF's filename?
		codes = self.parse_filename(cif_path)
		# CHALLENGE: ToBACCo has "organic nodes" as defined by ToBACCo that won't be picked up by MOFFLES
		# Combinining _on_ with linkers will have to wait.  In the meantime, just consider the pure metal cases.
		# Also skip B-containing sym_13_mc_12 and sym_16_mc_6, or Si-containing sym_4_on_14
		# sym_24_mc_13 will also be incompatible with our current decomposition scheme
		# For now, also sym_3_mc_0, sym_4_mc_1, sym_8_mc_7, sym_8_mc_8
		# TODO: Combine _on_ with the linker somehow
		# TODO: Consider challenges with combining certain metal+nitrogen nodes with carboxylates.
		# We might need SMIRKS transforms for those structures, too.
		if not any(False, [x in self.mof_db['nodes'] for x in codes['nodes']]):
			# Skip structures with tricky nodes (undefined for now)
			linkers = []
			for node in codes['nodes']:
				smiles = self.mof_db['nodes'][node]
				if smiles not in linkers:
					linkers.append(smiles)
			smiles = self.mof_db['linkers'][codes['linker']]
			if smiles != "None":
				linkers.append(smiles)
			linkers.sort()

			topology = codes['topology']
			if topology.startswith('test.'):
				topology = topology[5:]
			# Generate a reference MOFFLES based on SBU composition
			return assemble_moffles(linkers, topology, mof_name=codes['name'])
		else:
			return None

	def test_cif(self, cif_path):
		# Compares an arbitrary CIF file against its expected specification
		# Returns a formatted JSON string with the results
		# TODO: move to the parent class
		# TODO: consider updating compare_moffles with more helpful output data, including the name?
		moffles_from_name = self.expected_moffles(cif_path)
		if moffles_from_name is None:  # missing SBU info in the DB file
			return None
		else:
			# Calculate the MOFFLES derived from the CIF structure itself
			moffles_auto = cif2moffles(cif_path)
			return compare_moffles(moffles_from_name, moffles_auto, ['from_name', 'from_cif'])


if __name__ == "__main__":
	args = sys.argv[1:]
	if len(args) == 0:
		inputs = glob.glob(DEFAULT_CIFS + '/*.[Cc][Ii][Ff]')
	else:
		inputs = args

	# TODO: parse file type based on (future) command line argument, or guess it based on the filename
	# For now, just assume the file is a ToBACCo MOF
	comparer = TobaccoMOFs()
	moffles_results = []

	for cif_file in inputs:
		result = comparer.test_cif(cif_file)
		if result is not None:
			moffles_results.append(result)

	print summarize(moffles_results)
