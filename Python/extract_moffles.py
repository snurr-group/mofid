#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Calculate a preliminary MOFFLES code

Extract the node and linker identities from bin/sbu (my OpenBabel code) as
well as topological classification from Systre.  This wraps everything up
together.

@author: Ben Bucior
"""

import subprocess  # No support for timeouts except in Python 3.x, so do not use this for Systre
from easyprocess import EasyProcess  # Install with pip or from https://github.com/ponty/EasyProcess
import glob
import json
# import re
# import openbabel  # for visualization only, since my changes aren't backported to the python library
import sys, os

def path_to_resource(resource):
	# Get the path to resources, such as the MOF DB's or C++ code, without resorting to hardcoded paths
	python_path = os.path.dirname(__file__)
	return os.path.join(python_path, resource)

# Some default settings for my computer.  Adjust these based on your configuration:
SYSTRE_TIMEOUT = 30  # maximum time to allow Systre to run (seconds), since it hangs on certain CGD files
SBU_SYSTRE_PATH = "Test/topology.cgd"
GAVROG_LOC = path_to_resource("../Resources/External/Systre-1.2.0-beta2.jar")
SBU_BIN = path_to_resource("../bin/sbu")
if sys.platform == "win32":
	JAVA_LOC = "C:/Program Files/Java/jre1.8.0_102/bin/java"
elif sys.platform.startswith("linux"):
	JAVA_LOC = "java"
else:
	raise ValueError("Unknown platform.  Please specify file paths in Python/extract_moffles.py")


def extract_linkers(mof_path):
	# Extract MOF decomposition information using a C++ code based on OpenBabel
	cpp_run = EasyProcess([SBU_BIN, mof_path]).call()
	cpp_output = cpp_run.stdout
	sys.stderr.write(cpp_run.stderr)  # Re-forward sbu.cpp errors
	if cpp_run.return_code:  # EasyProcess uses threads, so you don't have to worry about the entire code crashing
		fragments = ["ERROR"]
	else:
		fragments = cpp_output.strip().split("\n")
		fragments = [x.strip() for x in fragments]  # clean up extra tabs, newlines, etc.

	cat = None
	if "simplified net(s)" in fragments[-1]:
		cat = fragments.pop()[6]  # "Found x simplified net(s)"
		cat = str(int(cat) - 1)
		if cat == "-1":
			cat = None

	return (sorted(fragments), cat)

def extract_topology(mof_path):
	# Extract underlying MOF topology using Systre and the output data from my C++ code
	java_run = EasyProcess([JAVA_LOC, "-Xmx1024m", "-cp", GAVROG_LOC, "org.gavrog.apps.systre.SystreCmdline", mof_path]).call(timeout=SYSTRE_TIMEOUT)
	java_output = java_run.stdout
	if java_run.timeout_happened:
		return "TIMEOUT"

	topologies = []  # What net(s) are found in the simplified framework(s)?
	current_component = 0
	topology_line = False
	for raw_line in java_output.split("\n"):
		line = raw_line.strip()
		if topology_line:
			topology_line = False
			rcsr = line.split()
			assert rcsr[0] == "Name:"
			topologies.append(rcsr[1])
		elif "ERROR" in line:
			return "ERROR"
		elif line == "Structure was identified with RCSR symbol:":
			topology_line = True
		elif line == "Structure is new for this run.":
			topologies.append("NEW")
		elif "Processing component " in line:
			assert len(topologies) == current_component  # Should extract one topology per component
			current_component += 1
			assert line[-2] == str(current_component)

	if len(topologies) == 0:
		return "ERROR"  # unexpected format
	first_net = topologies[0]  # Check that all present nets are consistent
	for net in topologies:
		if net != first_net:
			return "MISMATCH"
	return first_net

def assemble_moffles(linkers, topology, cat = None, mof_name="NAME_GOES_HERE"):
	# Assemble the MOFFLES code from its components
	moffles = ".".join(linkers) + " "
	moffles = moffles + "f1" + "."
	moffles = moffles + topology + "."
	if cat is not None:
		moffles = moffles + "cat" + cat + "."
	moffles = moffles + "F1" + "."
	moffles = moffles + mof_name
	return moffles

def parse_moffles(moffles):
	# Deconstruct a MOFFLES string into its pieces
	components = moffles.split()
	if len(components) == 1:
		if moffles.lstrip != moffles:  # Empty SMILES: no MOF found
			components.append(components[0])  # Move metadata to the right
			components[0] = ''
		else:
			raise ValueError("MOF metadata required")
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

def cif2moffles(cif_path):
	# Assemble the MOFFLES code from all of its pieces
	linkers, cat = extract_linkers(cif_path)
	if cat is not None:
		topology = extract_topology(SBU_SYSTRE_PATH)
	else:
		topology = "NA"
	mof_name = os.path.splitext(os.path.basename(cif_path))[0]
	return assemble_moffles(linkers, topology, cat, mof_name=mof_name)

def usage():
	raise SyntaxError("Run this code with a single parameter: path to the CIF file")

if __name__ == "__main__":
	args = sys.argv[1:]
	if len(args) != 1:
		usage()
	cif_file = args[0]
	
	print cif2moffles(cif_file)
