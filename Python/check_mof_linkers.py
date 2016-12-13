#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Calculate MOF linkers

Use the CSD criteria (no bonds to metals) in my modified OpenBabel code and
sbu.cpp to decompose MOFs into fragments.  Compare actual hMOF fragments
against the "recipe."

@author: Ben Bucior
"""

import subprocess
import glob
import json
# import re
import openbabel  # for visualization only, since my changes aren't backported to the python library

SBU_BIN = "C:/Users/Benjamin/Git/mofid/bin/sbu.exe"
CIF_DIR = "C:/Users/Benjamin/Desktop/Jiayi/Files/Dataset Comparison/hMOF"
HMOF_DB = "C:/Users/Benjamin/Git/mofid/Resources/hmof_linker_info.json"


def extract_linkers(mof_path):
	# Extract MOF decomposition information using a C++ code based on OpenBabel
	cpp_output = subprocess.check_output([SBU_BIN, mof_path])
	fragments = cpp_output.strip().split("\n")
	return sorted(fragments)

def parse_filename(hmof_path):
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

def load_hmof_components(db_file):
	# Load node and linker compositions from a saved JSON definition (to keep this file cleaner)
	with open(db_file, "r") as inp:
		hmof_db = json.load(inp)
	return hmof_db

def extract_db_smiles(db_dict):
	# Save all the nodes and linker definitions to the current directory
	# Then, convert to a big svg with
	# obabel linkers.can -O linkers.svg -xe -xC
	# Note: this function and the json file will be deprecated by the metal+linker
	# split paradigm suggested at group meeting ("molecule subtraction")
	linkers = db_dict["linkers"]
	for id in linkers:
		print(linkers[id] + " " + id)


if __name__ == "__main__":
	hmof_db = load_hmof_components(HMOF_DB)
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

