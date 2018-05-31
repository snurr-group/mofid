#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Removes metal nodes from a MOFFLES table

The names_to_tables script is powerful but still includes the metal nodes.
Let's remove those and generate yet another table, but this time only
including organic linkers (no metals).

@author: Ben Bucior
"""

import os, sys

import cheminformatics

def usage():
	print "Usage: python remove_metals.py smiles_part.tsv > smiles_nonmetals.tsv"
	exit()

# Establish a list of nonmetals based on sbu.cpp
# For now, let's disallow metals within the organic linkers, for simplicity
NONMETALS = [1, 2, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18, 32, 33, 34, 35, 36, 52, 53, 54, 85, 86]
def isMetal(atom_obj):
	# Operates like sbu.cpp:isMetal, but on a pybel atom instead of OBAtom*
	# The atom is classified as a metal if it's not a "nonmetal"
	return not(atom_obj.atomicnum in NONMETALS)

def contains_metal(smiles):
	# Does a SMILES entry contain a metal anywhere?
	mol = cheminformatics.pybel.readstring("smi", smiles)
	for atom in mol:
		if isMetal(atom):
			return True
	return False  # default if no metals found


if __name__ == "__main__":
	args = sys.argv[1:]
	if len(args) != 1:
		usage()
	
	with open(args[0], "r") as f:
		tsv_data = f.readlines()
		tsv_data = [x.rstrip("\n") for x in tsv_data]
	
	tsv_header = False
	if tsv_header:
		print tsv_data.pop(0)
	
	for line in tsv_data:
		smiles = line.split("\t")[1]
		if not (smiles=="*" or contains_metal(smiles)):
			print line
	
# On the output file, an interesting bash or SQL test you can run is
# the number of nonmetal linkers per MOF.  Are these two the same?
# wc -l OUTPUT/smiles_nonmetals.tsv
# cut -f1 OUTPUT/smiles_nonmetals.tsv | sort | uniq | wc -l
