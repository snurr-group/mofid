#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Calculate possible MOFFLES for a database

Use a list of CIF filenames to print the expected MOFFLES based on the
"recipe," mostly for the hMOFs or GA currently.

TODO: Consider extending this to other MOF databases, such as GA.
Though figuring out what to do with functionalization might be tricky.

@author: Ben Bucior
"""

import os, sys

from check_mof_linkers import TobaccoMOFs
from extract_moffles import assemble_moffles

def usage():
	print "Usage: python list_to_moffles.py FILENAME_LIST.txt > MOFFLES.smi"
	exit()

if __name__ == "__main__":
	args = sys.argv[1:]
	if len(args) != 1:
		usage()
	
	with open(args[0], "r") as f:
		file_names = f.readlines()
		file_names = [x.rstrip("\n") for x in file_names]
	
	compiler = TobaccoMOFs()  # TODO: make this more general to other DBs
	for line in file_names:
		expectation = compiler.expected_moffles(line)
		if expectation is None:
			print assemble_moffles("*", "NA", mof_name=compiler.parse_filename(line)['name'])
		else:
			print compiler.expected_moffles(line)['default']
	
