"""
Extracts metal composition based on MOFid

Based on remove_metals.py, but with the purpose of saving the metal compositions

@author: Ben Bucior
"""

import sys
import old_cheminformatics

def usage():
	print("Usage: python extract_metals.py smiles.tsv > smiles_metals.tsv")
	exit()

# Establish a list of nonmetals based on sbu.cpp
# For now, let's disallow metals within the organic linkers, for simplicity
NONMETALS = [1, 2, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18, 32, 33, 34, 35, 36, 52, 53, 54, 85, 86]
def isMetal(atom_obj):
	# Operates like sbu.cpp:isMetal, but on a pybel atom instead of OBAtom*
	# The atom is classified as a metal if it's not a "nonmetal"
	return not(atom_obj.atomicnum in NONMETALS)

def get_metals(smiles):
	# Get a list of metal names within a SMILES entry
	metals = []
	mol = old_cheminformatics.pybel.readstring("smi", smiles)
	for atom in mol:
		if isMetal(atom):
			metal_num = atom.atomicnum
			if metal_num not in metals:
				metals.append(metal_num)
	return metals

if __name__ == "__main__":
	args = sys.argv[1:]
	if len(args) != 1:
		usage()
	
	with open(args[0], "r") as f:
		tsv_data = f.readlines()
		tsv_data = [x.rstrip("\n") for x in tsv_data]
	
	tsv_header = False
	if tsv_header:
		print(tsv_data.pop(0))
	
	for line in tsv_data:
		[refcode, smiles] = line.split("\t")
		metals = get_metals(smiles)
		if len(metals) == 0:
			metals = ['NA']
		for metal in metals:
			print("\t".join([refcode, str(metal)]))
	
