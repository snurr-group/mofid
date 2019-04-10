'''
Parent module for obtaining MOFid data for a single .cif

@author: Ben Bucior
'''

import sys
import os
from mofid.id_constructor import (extract_fragments, extract_topology,
	assemble_mofkey, assemble_mofid, parse_mofid)
DEFAULT_OUTPUT_PATH = 'Output'

def cif2mofid(cif_path,output_path=DEFAULT_OUTPUT_PATH):
	# Assemble the MOFid string from all of its pieces.
	# Also export the MOFkey in an output dict for convenience.
	fragments, cat, base_mofkey = extract_fragments(cif_path,
		output_path)

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

	mofid = assemble_mofid(fragments, topology, cat, 
			mof_name=mof_name)
	parsed = parse_mofid(mofid)

	identifiers = {
		'mofid' : mofid,
		'mofkey' : mofkey,
		'smiles' : parsed['smiles'],
		'topology' : parsed['topology'],
		'cat' : parsed['cat'],
		'cifname' : parsed['name']
	}

	with open(os.path.join(output_path, 'python_mofid.txt'), 'w') as f:
		f.write(identifiers['mofid'] + '\n')
	with open(os.path.join(output_path, 'python_mofkey.txt'), 'w') as f:
		f.write(identifiers['mofkey'] + '\n')

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

	# Write MOFid and MOFkey output to files.
	with open(os.path.join(output_path, 'python_mofid.txt'), 'w') as f:
		f.write(identifiers['mofid'] + '\n')
	with open(os.path.join(output_path, 'python_mofkey.txt'), 'w') as f:
		f.write(identifiers['mofkey'] + '\n')