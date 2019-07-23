"""
Parent module for obtaining MOFid data for a single .cif

@author: Ben Bucior
"""

import sys
import os
from mofid.id_constructor import (extract_fragments, extract_topology,
	assemble_mofkey, assemble_mofid, parse_mofid)
from mofid.cpp_cheminformatics import openbabel_GetSpacedFormula
DEFAULT_OUTPUT_PATH = 'Output'

def cif2mofid(cif_path,output_path=DEFAULT_OUTPUT_PATH):
	# Assemble the MOFid string from all of its pieces.
	# Also export the MOFkey in an output dict for convenience.
	cif_path = os.path.abspath(cif_path)
	output_path = os.path.abspath(output_path)

	node_fragments, linker_fragments, cat, base_mofkey = extract_fragments(
		cif_path, output_path)
	if cat is not None:
		sn_topology = extract_topology(os.path.join(output_path,
			'SingleNode','topology.cgd'))
		an_topology = extract_topology(os.path.join(output_path,
			'AllNode','topology.cgd'))
		if sn_topology == an_topology:
			topology = sn_topology
		else:
			topology = sn_topology + ',' + an_topology
	else:
		topology = 'NA'

	mof_name = os.path.splitext(os.path.basename(cif_path))[0]
	mofkey = base_mofkey

	if topology != 'NA':
		base_topology = topology.split(',')[0]
		mofkey = assemble_mofkey(mofkey, base_topology)

	all_fragments = []
	all_fragments.extend(node_fragments)
	all_fragments.extend(linker_fragments)
	all_fragments.sort()
	mofid = assemble_mofid(all_fragments, topology, cat, 
			mof_name=mof_name)
	parsed = parse_mofid(mofid)

	identifiers = {
		'mofid' : mofid,
		'mofkey' : mofkey,
		'smiles_nodes' : node_fragments,
		'smiles_linkers' : linker_fragments,
		'smiles' : parsed['smiles'],
		'topology' : parsed['topology'],
		'cat' : parsed['cat'],
		'cifname' : parsed['name']
	}
	
	# Write MOFid and MOFkey output to files, as well as node/linker info
	with open(os.path.join(output_path, 'python_mofid.txt'), 'w') as f:
		f.write(identifiers['mofid'] + '\n')
	with open(os.path.join(output_path, 'python_mofkey.txt'), 'w') as f:
		f.write(identifiers['mofkey'] + '\n')
	with open(os.path.join(output_path, 'python_smiles_parts.txt'), 'w') as f:
		for smiles in node_fragments:
			f.write('node' + '\t' + smiles + '\n')
		for smiles in linker_fragments:
			f.write('linker' + '\t' + smiles + '\n')
	with open(os.path.join(output_path, 'python_molec_formula.txt'), 'w') as f:
		f.write(
			openbabel_GetSpacedFormula(
				os.path.join(output_path, 'orig_mol.cif'), ' ', False
			) + '\n'
		)

	return identifiers

if __name__ == '__main__':
	args = sys.argv[1:]
	if len(args) != 1 and len(args) != 2:
		raise SyntaxError('Usage: python run_mofid.py path_to_cif_for_analysis.cif OutputPathIfNonstandard')
	cif_file = args[0]
	output_path = DEFAULT_OUTPUT_PATH
	if len(args) == 2:
		output_path = args[1]

	identifiers = cif2mofid(cif_file, output_path)
	print(identifiers['mofid'])
	#print(identifiers['mofkey'])  # but incompatible with the use of stdout in run_folder.sh
