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
from extract_moffles import cif2moffles, assemble_moffles

SBU_BIN = "C:/Users/Benjamin/Git/mofid/bin/sbu.exe"
CIF_DIR = "C:/Users/Benjamin/Desktop/Jiayi/Files/Dataset Comparison/hMOF"
HMOF_DB = "C:/Users/Benjamin/Git/mofid/Resources/hmof_linker_info.json"
TOBACCO_DB = "C:/Users/Benjamin/Git/mofid/Resources/tobacco_info.json"
MOF_TYPE = "tobacco"  # hmof or tobacco


def any(member, list):
	# Is the member any part of the list?
	return member in list

def extract_linkers(mof_path):
	# Extract MOF decomposition information using a C++ code based on OpenBabel
	cpp_output = subprocess.check_output([SBU_BIN, mof_path])
	fragments = cpp_output.strip().split("\n")
	return sorted(fragments)

def parse_hmof_name(hmof_path):
	# Extract hMOF recipes from the filename, formatted as xxxhypotheticalMOF_####_i_#_j_#_k_#_m_#.cif
	codes = {"num": None, "i": None, "j": None, "k": None, "m": None, "cat": 0}
	parts = hmof_path.strip(".cifCIF").split("_")
	flag = None
	for part in parts:
		if flag is not None:
			#codes[flag] = int(part)
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

def parse_tobacco_name(tobacco_path):
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

	#print tobacco_path, codes
	return codes

def load_components(db_file):
	# Load node and linker compositions from a saved JSON definition (to keep this file cleaner)
	with open(db_file, "r") as inp:
		mof_db = json.load(inp)
	return mof_db

def extract_db_smiles(db_dict):
	# Save all the nodes and linker definitions to the current directory
	# Then, convert to a big svg with
	# obabel linkers.can -O linkers.svg -xe -xC
	# Note: this function and the json file will be deprecated by the metal+linker
	# split paradigm suggested at group meeting ("molecule subtraction")
	linkers = db_dict["linkers"]
	for id in linkers:
		print(linkers[id] + " " + id)

# Set filename parser
if MOF_TYPE == "hmof":
	parse_filename = parse_hmof_name
	SBU_DB = HMOF_DB
elif MOF_TYPE == "tobacco":
	parse_filename = parse_tobacco_name
	SBU_DB = TOBACCO_DB
else:
	raise ValueError("MOF parsing only (partially) implemented for hMOFs and ToBACCo")


if __name__ == False:  # disable my old hMOF code for now.  Re-incorporate later
	hmof_db = load_components(SBU_DB)
	# print hmof_db["nodes"]["0"][1]
	new_linkers = []
	for cif_file in glob.glob(CIF_DIR + '/*.[Cc][Ii][Ff]'):
		id = parse_filename(cif_file)
		if id['m'] != "0" or id["i"]!="0" or id["k"]!=id["j"] or id["j"] in new_linkers:
			# Later, this will probably just be id['m']!='0' since we'll want to temporarily exclude functionalization
			continue
		print cif_file
		for fragment in extract_linkers(cif_file):
			print fragment
		new_linkers.append(id["j"])
		print "\n"
	extract_db_smiles(hmof_db)



if __name__ == "__main__":
	mof_db = load_components(SBU_DB)
	for cif_file in glob.glob('C:/Users/Benjamin/Desktop/ToBACCo - Copy/output_structures' + '/*.[Cc][Ii][Ff]'):
		id = parse_filename(cif_file)
		if id['nodes'] == ['sym_6_mc_3'] and id['linker'] == "L_12" and id['topology'] == "test.pcu":
			print "Found MOF-5!:", id
		# CHALLENGE: ToBACCo has "organic nodes" as defined by ToBACCo that won't be picked up by MOFFLES
		# Combinining _on_ with linkers will have to wait.  In the meantime, just consider the pure metal cases.
		# Also skip B-containing sym_13_mc_12 and sym_16_mc_6, or Si-containing sym_4_on_14
		# sym_24_mc_13 will also be incompatible with our current decomposition scheme
		# For now, also sym_3_mc_0, sym_4_mc_1, sym_8_mc_7, sym_8_mc_8
		# TODO: Combine _on_ with the linker somehow
		#if not any(True, ['_on_' in x for x in id['nodes']]):
		if not any(False, [x in mof_db['nodes'] for x in id['nodes']]):
			#print id
			linkers = []
			for node in id['nodes']:
				smiles = mof_db['nodes'][node]
				if smiles not in linkers:
					linkers.append(smiles)
			smiles = mof_db['linkers'][id['linker']]
			if smiles != "None":
				linkers.append(smiles)
			linkers.sort()
			#print linkers

			topology = id['topology']
			if topology.startswith('test.'):
				topology = topology[5:]

			moffles_name = assemble_moffles(linkers, topology, mof_name=id['name'])
			# print moffles_name
			#moffles_auto = cif2moffles(cif_file.replace("\\", "/"))
			moffles_auto = cif2moffles(cif_file)
			if moffles_auto == moffles_name:
				print "Success!:", moffles_auto
			else:
				print "Failure::"
				print "Database:", moffles_name
				print "Auto:", moffles_auto
			# Need to call the extract_moffles code #cif2moffles
			# Generate a comparison MOFFLES based on SBU composition

	
