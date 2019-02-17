"""
Calculate possible MOFFLES for a database

Use a list of CIF filenames to print the expected MOFFLES based on the
"recipe," mostly for the hMOFs or GA currently.

Long-term, we could consider extending this to other MOF databases, such as
the GA MOFs.  But there are multiple things that could prove difficult about
the implementation: functionalization is inexact, the notation/positioning
for nitrogen-terminated vs. carboxylate linkers is indeterminate, topologies
aren't necessarily what you would expect, etc.  So let's focus on the ToBaCCo
MOFs for generating linker SMILES, for now.

@author: Ben Bucior
"""

import sys
from check_mof_linkers import TobaccoMOFs
from run_mofid import assemble_moffles

if __name__ == "__main__":
	args = sys.argv[1:]
	if len(args) != 1:
		raise SyntaxError('Usage: python list_to_mofid.py FILENAME_LIST.txt > MOFFLES.smi')
	
	with open(args[0], "r") as f:
		file_names = f.readlines()
		file_names = [x.rstrip("\n") for x in file_names]
	
	compiler = TobaccoMOFs()  # TODO: make this more general to other DBs
	for line in file_names:
		expectation = compiler.expected_moffles(line)
		if expectation is None:
			print(assemble_moffles("*", "NA",
				mof_name=compiler.parse_filename(line)['name']))
		else:
			print(compiler.expected_moffles(line)['default'])
	
