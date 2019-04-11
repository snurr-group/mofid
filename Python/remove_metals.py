"""
Removes metal nodes from a MOFid table

The names_to_tables script is powerful but still includes the metal nodes.
Let's remove those and generate yet another table, but this time only
including organic linkers (no metals).

@author: Ben Bucior
"""

import sys
import mofid.extract_metals as extract_metals

def contains_metal(smiles):
	# Does a SMILES entry contain a metal anywhere?
	if len(extract_metals.get_metals(smiles)) > 0:
		return True
	else:
		return False


if __name__ == '__main__':
	args = sys.argv[1:]
	if len(args) != 1:
		raise SyntaxError('Usage: python remove_metals.py smiles_part.tsv > smiles_nonmetals.tsv')
	
	with open(args[0], 'r') as f:
		tsv_data = f.readlines()
		tsv_data = [x.rstrip('\n') for x in tsv_data]
	
	tsv_header = False
	if tsv_header:
		print(tsv_data.pop(0))
	
	for line in tsv_data:
		smiles = line.split('\t')[1]
		if not (smiles=='*' or contains_metal(smiles)):
			print(line)
	
# On the output file, an interesting bash or SQL test you can run is
# the number of nonmetal linkers per MOF.  Are these two the same?
# wc -l OUTPUT/smiles_nonmetals.tsv
# cut -f1 OUTPUT/smiles_nonmetals.tsv | sort | uniq | wc -l
