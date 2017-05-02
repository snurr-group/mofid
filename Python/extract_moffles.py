#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Calculate a preliminary MOFFLES code

Extract the node and linker identities from bin/sbu (my OpenBabel code) as
well as topological classification from Systre.  This wraps everything up
together.

@author: Ben Bucior
"""

import subprocess
import glob
import json
# import re
# import openbabel  # for visualization only, since my changes aren't backported to the python library
import sys, os

# Some default settings for my computer.  Adjust these based on your configuration:
SBU_SYSTRE_PATH = "Test/topology.cgd"
if sys.platform == "win32":
	SBU_BIN = "C:/Users/Benjamin/Git/mofid/bin/sbu.exe"
	# Settings for Systre
	JAVA_LOC = "C:/Program Files/Java/jre1.8.0_102/bin/java"
	GAVROG_LOC = "C:/Users/Benjamin/Software/Gavrog-0.6.0/Systre.jar"
elif sys.platform.startswith("linux"):
	SBU_BIN = "/home/bbucior/Git/mofid/bin/sbu"
	JAVA_LOC = "java"
	GAVROG_LOC = "/home/bbucior/Software/Gavrog-0.6.0/Systre.jar"
else:
	raise ValueError("Unknown platform.  Please specify file paths in Python/extract_moffles.py")


def extract_linkers(mof_path):
	# Extract MOF decomposition information using a C++ code based on OpenBabel
	cpp_output = subprocess.check_output([SBU_BIN, mof_path])
	fragments = cpp_output.strip().split("\n")
	fragments = [x.strip() for x in fragments]  # clean up extra tabs, newlines, etc.
	return sorted(fragments)

def extract_topology(mof_path):
	# Extract underlying MOF topology using Systre and the output data from my C++ code
	java_output = subprocess.check_output([JAVA_LOC, "-Xmx512m", "-cp", GAVROG_LOC, "org.gavrog.apps.systre.SystreCmdline", mof_path])
	# Put the BATCH ONE HERE
	topology_line = False
	for line in java_output.split("\n"):
		if topology_line:
			topology_line = False
			rcsr = line.strip().split()
			assert rcsr[0] == "Name:"
			return rcsr[1]
		elif line.strip() == "Structure was identified with RCSR symbol:":
			topology_line = True
		elif line.strip() == "Structure is new for this run.":
			return "NEW"
	return "ERROR"  # unexpected format

def assemble_moffles(linkers, topology, cat = "CAT_TBD", mof_name="NAME_GOES_HERE"):
	# Assemble the MOFFLES code from its components
	moffles = ".".join(linkers) + " "
	moffles = moffles + "f1" + "."
	moffles = moffles + topology + "."
	moffles = moffles + cat + "."
	moffles = moffles + "F1" + "."
	moffles = moffles + mof_name
	return moffles

def parse_moffles(moffles):
	# Deconstruct a MOFFLES string into its pieces
	components = moffles.split()
	smiles = components[0]
	if len(components) > 2:
		print "Bad MOFFLES:", moffles
		raise ValueError("FIXME: parse_moffles currently does not support spaces in common names")
	metadata = components[1]
	metadata = metadata.split('.')

	mof_name = None
	cat = None
	topology = None
	for loc, tag in enumerate(metadata):
		if loc == 0 and tag != 'f1':
			raise ValueError("MOFFLES v1 must start with tag 'f1'")
		if tag == 'F1':
			mof_name = '.'.join(metadata[loc+1:])
			break
		elif tag.lower().startswith('cat'):
			cat = tag[3:]
		else:
			topology = tag
	return dict(
		smiles = smiles,
		topology = topology,
		cat = cat,
		name = mof_name
	)

def compare_moffles(moffles1, moffles2, names=None):
	# Compares MOFFLES strings to identify sources of difference, if any
	if names is None:
		names = ['mof1', 'mof2']
	if moffles1 is None or moffles2 is None:
		mof_name = 'Undefined'
		for x in [moffles1, moffles2]:
			if x is not None:
				mof_name = parse_moffles(x)['name']
		return {'match': 'NA',
		        'errors': ['Undefined composition'],
		        'topology': None,
		        'smiles': None,
		        'cat': None,
		        names[0]: moffles1,
		        names[1]: moffles2,
		        'name': mof_name
		        }
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
			comparison['errors'].append("err_" + key)
	return comparison

def cif2moffles(cif_path):
	# Assemble the MOFFLES code from all of its pieces
	linkers = extract_linkers(cif_path)
	topology = extract_topology(SBU_SYSTRE_PATH)
	mof_name = os.path.splitext(os.path.basename(cif_path))[0]
	return assemble_moffles(linkers, topology, mof_name=mof_name)

def usage():
	raise SyntaxError("Run this code with a single parameter: path to the CIF file")

if __name__ == "__main__":
	args = sys.argv[1:]
	if len(args) != 1:
		usage()
	cif_file = args[0]
	
	print cif2moffles(cif_file)
