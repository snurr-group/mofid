"""
Calculate a MOFid string

Extract the node and linker identities from bin/sbu (my Open Babel code) as
well as topological classification from Systre.  This script runs everything
to wrap it up as a formated string.

@author: Ben Bucior
"""

import sys
import os
from mofid.paths import resources_path, bin_path

if sys.version_info[0] < 3:
	try:
		import subprocess32 as subprocess
	except:
		raise AssertionError('You must install subprocess32 if running Python2')
else:
	import subprocess

# Make sure Java is in user's path
try:
	subprocess.call('java',stderr=subprocess.PIPE)
except EnvironmentError:
	raise AssertionError('You must have Java in your path!')

# Make sure Java is in user's path
try:
	subprocess.call('java',stderr=subprocess.PIPE)
except EnvironmentError:
	raise AssertionError('You must have Java in your path!')

# Some default paths
GAVROG_LOC = os.path.join(resources_path,'Systre-1.2.0-beta2.jar')
JAVA_LOC = 'java'
RCSR_PATH = os.path.join(resources_path,'RCSRnets.arc')
DEFAULT_SYSTRE_CGD = os.path.join('Output','SingleNode','topology.cgd')
SYSTRE_TIMEOUT = 30  # max time to allow Systre to run (seconds), since it hangs on certain CGD files
SBU_BIN = os.path.join(bin_path,'sbu')

# Can update the RCSR version at http://rcsr.anu.edu.au/systre
SYSTRE_CMD_LIST = [
	JAVA_LOC,
	'-Xmx1024m',  # allocate up to 1GB of memory
	'-cp', GAVROG_LOC,  # call a specific classpath in the .jar file
	'org.gavrog.apps.systre.SystreCmdline',
	RCSR_PATH  # RCSR archive to supplement the old version in the .jar file
	]

def runcmd(cmd_list, timeout=None):
	if timeout is None:
		return subprocess.run(cmd_list, universal_newlines=True,
			stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	else:
		return subprocess.run(cmd_list, universal_newlines=True,
			stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=timeout)

def extract_fragments(mof_path,output_path):
	# Extract MOF decomposition information using a C++ code based on OpenBabel
	cpp_run = runcmd([SBU_BIN, mof_path, output_path])
	cpp_output = cpp_run.stdout
	sys.stderr.write(cpp_run.stderr)  # Re-forward sbu.cpp errors
	if cpp_run.returncode:  # EasyProcess uses threads, so you don't have to worry about the entire code crashing
		fragments = ['*']  # Null-behaving atom for Open Babel and rdkit, so the .smi file is still useful
	else:
		all_fragments = cpp_output.strip().split('\n')
		all_fragments = [x.strip() for x in all_fragments]  # clean up extra tabs, newlines, etc.

	cat = None
	if 'simplified net(s)' in all_fragments[-1]:
		cat = all_fragments.pop()[8]  # '# Found x simplified net(s)'
		cat = str(int(cat) - 1)
		if cat == '-1':
			cat = None
	
	# Parse node/linker fragment notation
	if all_fragments[0] != '# Nodes:':
		return (None, None, cat, None)
	all_fragments.pop(0)
	linker_flag_loc = all_fragments.index('# Linkers:')
	node_fragments = all_fragments[:linker_flag_loc]
	linker_fragments = all_fragments[linker_flag_loc+1:]  # could be the empty set

	base_mofkey = None
	if not cpp_run.returncode:  # If it's a successful run
		mofkey_loc = os.path.join(
			output_path, 'MetalOxo', 'mofkey_no_topology.txt')
		with open(mofkey_loc) as f:
			base_mofkey = f.read().rstrip()  # ending newlines, etc.

	return (sorted(node_fragments), sorted(linker_fragments), cat, base_mofkey)

def extract_topology(mof_path):
	# Extract underlying MOF topology using Systre and the output data from my C++ code
	try:
		java_run = runcmd(SYSTRE_CMD_LIST + [mof_path],
			timeout=SYSTRE_TIMEOUT)
	except subprocess.TimeoutExpired:
		return 'TIMEOUT'
	java_output = java_run.stdout

	topologies = []  # What net(s) are found in the simplified framework(s)?
	current_component = 0
	topology_line = False
	repeat_line = False
	for raw_line in java_output.split('\n'):
		line = raw_line.strip()
		if topology_line:
			topology_line = False
			rcsr = line.split()
			assert rcsr[0] == 'Name:'
			topologies.append(rcsr[1])
		elif repeat_line:
			repeat_line = False
			assert line.split()[0] == 'Name:'
			components = line.split('_')  # Line takes the form 'Name:    refcode_clean_component_x'
			assert components[-2] == 'component'
			topologies.append(topologies[int(components[-1]) - 1])  # Subtract one since Systre is one-indexed
		elif 'ERROR' in line:
			return 'ERROR'
		elif 'Structure was found in archive' in line:
			topology_line = True
		elif line == 'Structure is new for this run.':
			# This line is only printed if new to both versions of the RCSR database:
			# a copy saved in the .jar file and an updated version in Resources/RCSRnets.arc
			topologies.append('UNKNOWN')
		elif line == 'Structure already seen in this run.':
			repeat_line = True
		elif 'Processing component ' in line:
			assert len(topologies) == current_component  # Should extract one topology per component
			current_component += 1
			assert line[-2] == str(current_component)

	if len(topologies) == 0:
		return 'ERROR'  # unexpected format
	first_net = topologies[0]  # Check that all present nets are consistent
	for net in topologies:
		if net != first_net:
			return 'MISMATCH'
	return first_net

def assemble_mofid(fragments, topology, cat = None, mof_name='NAME_GOES_HERE'):
	# Assemble the MOFid string from its components
	mofid = '.'.join(fragments) + ' '
	mofid = mofid + 'MOFid-v1' + '.'
	mofid = mofid + topology + '.'
	if cat == 'no_mof':
		mofid = mofid + cat
	elif cat is not None:
		mofid = mofid + 'cat' + cat
	if mofid.startswith(' '):  # Null linkers.  Make .smi compatible
		mofid = '*' + mofid + 'no_mof'
	mofid = mofid + ';' + mof_name
	return mofid

def assemble_mofkey(base_mofkey, base_topology):
	# Add a topology to an existing MOFkey
	return base_mofkey.replace('MOFkey-v1', 'MOFkey-v1.' + base_topology)

def parse_mofid(mofid):
	# Deconstruct a MOFid string into its pieces
	#[mofid_data, mofid_name]
	# Normalize the string by removing trailing spaces, e.g. newlines
	mofid_parts = mofid.rstrip().split(';')
	mofid_data = mofid_parts[0]
	if len(mofid_parts) > 1:
		mofid_name = ';'.join(mofid_parts[1:])
	else:
		mofid_name = None

	components = mofid_data.split()
	if len(components) == 1:
		if mofid_data.lstrip != mofid_data:  # Empty SMILES: no MOF found
			components.append(components[0])  # Move metadata to the right
			components[0] = ''
		else:
			raise ValueError('MOF metadata required')
	smiles = components[0]
	if len(components) > 2:
		raise ValueError('Bad MOFid containing extra spaces before the semicolon:' + mofid)
	metadata = components[1]
	metadata = metadata.split('.')

	cat = None
	topology = None
	for loc, tag in enumerate(metadata):
		if loc == 0 and not tag.startswith('MOFid'):
			raise ValueError('MOFid-v1 must start with the correct tag')
		if loc == 0 and tag[5:] != '-v1':
			raise ValueError('Unsupported version of MOFid')
		elif loc == 1:
			topology = tag
		elif tag.lower().startswith('cat'):
			cat = tag[3:]
		else:
			pass  # ignoring other MOFid tags, at least for now
	return dict(
		smiles = smiles,
		topology = topology,
		cat = cat,
		name = mofid_name
	)
