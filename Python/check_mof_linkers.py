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
# import re
import openbabel  # for visualization only, since my changes aren't backported to the python library
from extract_moffles import cif2moffles, assemble_moffles, parse_moffles

SBU_BIN = "C:/Users/Benjamin/Git/mofid/bin/sbu.exe"
CIF_DIR = "C:/Users/Benjamin/Desktop/Jiayi/Files/Dataset Comparison/hMOF"
HMOF_DB = "C:/Users/Benjamin/Git/mofid/Resources/hmof_linker_info.json"
TOBACCO_DB = "C:/Users/Benjamin/Git/mofid/Resources/tobacco_info.json"
KNOWN_DB = "C:/Users/Benjamin/Git/mofid/Resources/known_info.json"
MOF_TYPE = "tobacco"  # hmof or tobacco


def any(member, list):
	# Is the member any part of the list?
	return member in list

class MOFCompare:
	# Compares a MOF's decomposition against its expected components
	def __init__(self):
		self.db_file = None
		self.load_components()
		raise UserWarning("MOF parsing only (partially) implemented for hMOFs and ToBACCo")

	def extract_cif_linkers(self, mof_path):
		# Extract MOF decomposition information using a C++ code based on OpenBabel
		cpp_output = subprocess.check_output([SBU_BIN, mof_path])
		fragments = cpp_output.strip().split("\n")
		return sorted(fragments)

	def compare_moffles(self, moffles1, moffles2, names=None):
		# Compares MOFFLES strings to identify sources of difference, if any
		if names is None:
			names = ['mof1', 'mof2']
		if moffles1 is None or moffles2 is None:
			# TODO: add MOF names, other fields?
			return {'match': 'NA', 'errors': ['Undefined composition'], 'topology': None, 'linkers': None}
		parsed = [parse_moffles(x) for x in [moffles1, moffles2]]
		comparison = dict()
		comparison['match'] = True
		comparison['errors'] = []
		comparison[names[0]] = moffles1
		comparison[names[1]] = moffles2
		for key in parsed[0]:
			expected = parsed[0][key]
			if parsed[1][key] == expected:
				comparison[key] = expected
			else:
				comparison[key] = False
				comparison['match'] = False
				comparison['errors'].append(key)
		return comparison

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


class HypoMOFs(MOFCompare):
	def __init__(self):
		self.db_file = HMOF_DB
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
		hmof_db = load_components(SBU_DB)
		# print hmof_db["nodes"]["0"][1]
		new_linkers = []
		for cif_file in glob.glob(CIF_DIR + '/*.[Cc][Ii][Ff]'):
			id = parse_filename(cif_file)
			if id['m'] != "0" or id["i"] != "0" or id["k"] != id["j"] or id["j"] in new_linkers:
				# Later, this will probably just be id['m']!='0' since we'll want to temporarily exclude functionalization
				continue
			print cif_file
			for fragment in extract_cif_linkers(cif_file):
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

	def test_cif(self, cif_path):
		# Compares an arbitrary CIF file against its expected specification
		# Returns a formatted JSON string with the results
		##### TODO: clean this up significantly
		id = self.parse_filename(cif_path)
		# CHALLENGE: ToBACCo has "organic nodes" as defined by ToBACCo that won't be picked up by MOFFLES
		# Combinining _on_ with linkers will have to wait.  In the meantime, just consider the pure metal cases.
		# Also skip B-containing sym_13_mc_12 and sym_16_mc_6, or Si-containing sym_4_on_14
		# sym_24_mc_13 will also be incompatible with our current decomposition scheme
		# For now, also sym_3_mc_0, sym_4_mc_1, sym_8_mc_7, sym_8_mc_8
		# TODO: Combine _on_ with the linker somehow
		# TODO: encapsulate code below into a relevant function
		moffles_name = None
		if not any(False, [x in self.mof_db['nodes'] for x in id['nodes']]):
			linkers = []
			for node in id['nodes']:
				smiles = self.mof_db['nodes'][node]
				if smiles not in linkers:
					linkers.append(smiles)
			smiles = self.mof_db['linkers'][id['linker']]
			if smiles != "None":
				linkers.append(smiles)
			linkers.sort()

			topology = id['topology']
			if topology.startswith('test.'):
				topology = topology[5:]
			# Generate a reference MOFFLES based on SBU composition
			moffles_name = assemble_moffles(linkers, topology, mof_name=id['name'])
			# TODO: Below: used to be indented one level out, but that required evaluating unknown structures
			# TODO: check for null comparison.  Maybe make a simple NULL compare_moffles that doesn't require moffles_auto eval?
			# Calculate the MOFFLES derived from the CIF structure itself
			moffles_auto = cif2moffles(cif_path)
			return self.compare_moffles(moffles_name, moffles_auto, ['from_name', 'from_cif'])
		return None


if __name__ == "__main__":
	args = sys.argv[1:]
	if len(args) == 0:
		inputs = glob.glob('Data/tobacco_L_12/quick' + '/*.[Cc][Ii][Ff]')
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

	# TODO: encapsulate this code in the relevant
	print moffles_results
	if True:  # Errors analysis
		print "Error analysis:"
		error_types = {'topology': 0, 'smiles': 0, 'both': 0, 'success': 0}
		for match in moffles_results:
			if match['match']:
				error_types['success'] += 1
			if len(match['errors']) > 1:
				error_types['both'] +=1
			elif 'topology' in match['errors']:
				error_types['topology'] += 1
			elif 'linkers' in match['errors']:
				error_types['smiles'] += 1
		# TODO: Handle the NA (error of "Undefined composition")
		print error_types, "Total:", len(moffles_results)
