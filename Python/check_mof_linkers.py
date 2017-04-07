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
import time
# import copy  # copy.deepcopy(x)
# import re
import openbabel  # for visualization only, since my changes aren't backported to the python library
from extract_moffles import cif2moffles, assemble_moffles, parse_moffles, compare_moffles, extract_linkers

# Locations of important files, relative to the Python source code
SBU_BIN = "../bin/sbu.exe"

HMOF_DB = "../Resources/hmof_linker_info.json"
TOBACCO_DB = "../Resources/tobacco_info.json"
KNOWN_DB = "../Resources/known_mof_info.json"

KNOWN_DEFAULT_CIFS = "../Data/RingCIFs"
TOBACCO_DEFAULT_CIFS = "../Data/tobacco_L_12/quick"
HMOF_DEFAULT_CIFS = "../Data/hmofs_i_0_no_cat"
NO_ARG_CIFS = KNOWN_DEFAULT_CIFS  # KnownMOFs() comparisons are used if no args are specified.  See arg parsing of main
PRINT_CURRENT_MOF = True


def any(member, list):
	# Is the member any part of the list?
	return member in list

def path_to_resource(resource):
	# Get the path to resources, such as the MOF DB's or C++ code, without resorting to hardcoded paths
	# However, some absolute paths are still present in extract_moffles.py since they're system-wide
	python_path = os.path.dirname(__file__)
	return os.path.join(python_path, resource)

def basename(path):
	# Get the basename for a given path, without the file extension
	return os.path.splitext(os.path.basename(path))[0]

def mof_log(msg):
	# Logging helper function, which writes to stderr to avoid polluting the json
	if PRINT_CURRENT_MOF:
		sys.stderr.write(msg)

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
	summarized['errors']['elapsed_time'] = sum([mof['time'] for mof in results])
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

	def test_cif(self, cif_path):
		# Compares an arbitrary CIF file against its expected specification
		# Returns a formatted JSON string with the results
		start = time.time()
		moffles_from_name = self.expected_moffles(cif_path)
		if moffles_from_name is None:  # missing SBU info in the DB file
			return None  # Currently, skip reporting of structures with undefined nodes/linkers
		else:
			# Calculate the MOFFLES derived from the CIF structure itself
			moffles_auto = cif2moffles(cif_path)
			comparison = compare_moffles(moffles_from_name, moffles_auto, ['from_name', 'from_cif'])
			comparison['time'] = time.time() - start
			comparison['name_parser'] = self.__class__.__name__
			return comparison


class KnownMOFs(MOFCompare):
	# Minimal class which doesn't have to do much work to scour the database of known MOFs.
	# Excellent as a integration test for my code, i.e. did my changes cause anything else to obviously break?
	def __init__(self):
		self.db_file = path_to_resource(KNOWN_DB)
		self.load_components()

	def parse_filename(self, mof_path):
		# Extract basename of the MOF.  expected_moffles will convert it to the reference moffles string
		return basename(mof_path)

	def expected_moffles(self, cif_path):
		# What is the expected MOFFLES based on the information in a MOF's filename?
		mof_name = self.parse_filename(cif_path)
		if mof_name in self.mof_db:
			return self.mof_db[mof_name]
		else:
			return None


class HypoMOFs(MOFCompare):
	def __init__(self):
		self.db_file = path_to_resource(HMOF_DB)
		self.load_components()

	def parse_filename(self, hmof_path):
		# Extract hMOF recipes from the filename, formatted as xxxhypotheticalMOF_####_i_#_j_#_k_#_m_#.cif
		codes = {"num": None, "i": None, "j": None, "k": None, "m": None, "cat": 0}
		mof_name = basename(hmof_path)  # Get the basename without file extension
		parts = mof_name.split("_")
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
		codes["name"] = mof_name
		return codes

	def expected_moffles(self, cif_path):
		# What is the expected MOFFLES based on the information in a MOF's filename?
		# TODO: add catentation and expected (known?) topology
		# For now, just assume everything is **pcu** since we're starting with Zn4O nodes
		codes = self.parse_filename(cif_path)
		code_key = {
			'i': 'nodes',
			'j': 'linkers',
			'k': 'linkers',
			'm': 'functionalization',
		}

		if codes['m'] != "0" or codes["k"] != codes["j"]:
			return None  # For now, temporarily exclude functionalization (until we can figure out SMIRKS transforms)
			# Also exclude multiple linkers, to narrow down the number of test structures, at least to begin

		is_component_defined = []
		for key in code_key:
			is_component_defined.extend([x in self.mof_db[code_key[key]] for x in codes[key]])

		if not any(False, is_component_defined):
			sbus = []
			sbu_codes = ['i', 'j', 'k']
			for part in sbu_codes:
				smiles = self.mof_db[code_key[part]][codes[part]]
				# print "SMILES:", smiles
				if smiles not in sbus:
					sbus.append(smiles)
			sbus.sort()

			# topology = "pcu"  # FIXME: temporary assumption for Zn4O nodes
			topology = "ERROR"  # Temporarily disable topology checks for hMOFs, since pcu is still buggy
			return assemble_moffles(sbus, topology, mof_name=codes['name'])
		else:
			return None


class TobaccoMOFs(MOFCompare):
	def __init__(self):
		self.db_file = path_to_resource(TOBACCO_DB)
		self.load_components()

	def parse_filename(self, tobacco_path):
		# Extract ToBACCo recipes from the filename
		# Format: topology_sym_x_node_type_sym_x_node2_type_L_linkernum.cif ("_" for empty linker)
		# Can we parse this using ToBACCo's own code??
		codes = {"name": None, "nodes": [], "linker": None, "topology": None}
		mof_info = basename(tobacco_path)  # Get the basename without file extension
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

class AutoCompare:
	# Automatically selects the appropriate MOF comparison class using filename-based heuristics
	# Note: this is not a full implementation of MOFCompare--just a wrapper for test_cif
	# TODO: allow the user to manually specify the MOF class using a command line flag
	def __init__(self):
		self.known = KnownMOFs()
		self.precalculated = self.known.mof_db.keys()
		self.hmof = HypoMOFs()
		self.tobacco = TobaccoMOFs()
		# Maybe also a NullMOFs class eventually, which is just a calculator sans comparisons?

	def test_cif(self, cif_path):
		# Dispatch to the class corresponding to the source of the input CIF
		mof_info = basename(cif_path)
		if mof_info in self.precalculated:
			mof_log("...using precompiled table of known MOFs\n")
			return self.known.test_cif(cif_path)
		elif "hypotheticalmof" in mof_info.lower() or "hmof" in mof_info.lower():
			mof_log("...parsing file with rules for hypothetical MOFs\n")
			return self.hmof.test_cif(cif_path)
		elif "_sym_" in mof_info:
			mof_log("...parsing file with rules for ToBACCo MOFs\n")
			return self.tobacco.test_cif(cif_path)
		else:
			mof_log("...unable to find a suitable rule automatically\n")
			return None


if __name__ == "__main__":
	comparer = AutoCompare()  # By default, guess the MOF type by filename
	args = sys.argv[1:]
	if len(args) == 0:  # validation testing against reference MOFs
		inputs = glob.glob(path_to_resource(NO_ARG_CIFS) + '/*.[Cc][Ii][Ff]')
		comparer = KnownMOFs()
	elif len(args) == 1 and args[0].endswith("/"):
		# Run a whole directory if specified as a single argument with an ending slash
		inputs = glob.glob(args[0] + '*.[Cc][Ii][Ff]')
	else:
		inputs = args

	moffles_results = []
	for num_cif, cif_file in enumerate(inputs):
		mof_log(" ".join(["Found CIF", str(num_cif+1), "of", str(len(inputs)), ":", cif_file]) + "\n")
		result = comparer.test_cif(cif_file)
		if result is not None:
			moffles_results.append(result)

	json.dump(summarize(moffles_results), sys.stdout, indent=4)
