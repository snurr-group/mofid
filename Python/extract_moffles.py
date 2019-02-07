"""
Calculate a preliminary MOFFLES code

Extract the node and linker identities from bin/sbu (my OpenBabel code) as
well as topological classification from Systre.  This wraps everything up
together.

@author: Ben Bucior
"""

import sys
import os

# Ensure that subprocess32 is installed if running Py2
if sys.version_info[0] < 3:
	try:
		import subprocess32 as subprocess
	except:
		raise AssertionError('You must install subprocess if running Python2')
else:
	import subprocess
	
# Make sure Java is in user's path
try:
	subprocess.call('java',stderr=subprocess.PIPE)
except EnvironmentError:
	raise AssertionError('You must have Java in your path!')

def path_to_resource(resource):
	# Get the path to resources, such as the MOF DB's or C++ code
	# without resorting to hardcoded paths
	python_path = os.path.dirname(__file__)
	return os.path.join(python_path, resource)

def runcmd(cmd_list, timeout=None):
	if timeout is None:
		return subprocess.run(cmd_list, universal_newlines=True,
			stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	else:
		return subprocess.run(cmd_list, universal_newlines=True,
			stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=timeout)


# Some default settings for my computer.  Adjust these based on your config:
SYSTRE_TIMEOUT = 30  # max time to allow Systre to run (seconds), since it hangs on certain CGD files
GAVROG_LOC = path_to_resource("../Resources/External/Systre-1.2.0-beta2.jar")
SBU_BIN = path_to_resource("../bin/sbu")
JAVA_LOC = "java"

# The default path is set in deconstructor.h:DEFAULT_OUTPUT_PATH.
# Update it here if the directory changes.
DEFAULT_OUTPUT_PATH = "Output/"
# Default to the single-node decomposition algorithm for assigning topology
DEFAULT_SYSTRE_CGD = os.path.join(DEFAULT_OUTPUT_PATH, "SingleNode/topology.cgd")


def extract_linkers(mof_path, output_file_path=DEFAULT_OUTPUT_PATH):
	# Extract MOF decomposition information using a C++ code based on OpenBabel
	cpp_run = runcmd([SBU_BIN, mof_path, output_file_path])
	cpp_output = cpp_run.stdout
	sys.stderr.write(cpp_run.stderr)  # Re-forward sbu.cpp errors
	if cpp_run.returncode:  # EasyProcess uses threads, so you don't have to worry about the entire code crashing
		fragments = ["*"]  # Null-behaving atom for Open Babel and rdkit, so the .smi file is still useful
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
	try:
		java_run = runcmd([JAVA_LOC, "-Xmx1024m", "-cp", GAVROG_LOC,
			"org.gavrog.apps.systre.SystreCmdline", mof_path],
			timeout=SYSTRE_TIMEOUT)
	except subprocess.TimeoutExpired:
		return "TIMEOUT"
	java_output = java_run.stdout

	topologies = []  # What net(s) are found in the simplified framework(s)?
	current_component = 0
	topology_line = False
	repeat_line = False
	for raw_line in java_output.split("\n"):
		line = raw_line.strip()
		if topology_line:
			topology_line = False
			rcsr = line.split()
			assert rcsr[0] == "Name:"
			topologies.append(rcsr[1])
		elif repeat_line:
			repeat_line = False
			assert line.split()[0] == "Name:"
			components = line.split("_")  # Line takes the form "Name:    refcode_clean_component_x"
			assert components[-2] == "component"
			topologies.append(topologies[int(components[-1]) - 1])  # Subtract one since Systre is one-indexed
		elif "ERROR" in line:
			return "ERROR"
		elif line == "Structure was identified with RCSR symbol:":
			topology_line = True
		elif line == "Structure is new for this run.":
			topologies.append("NEW")
		elif line == "Structure already seen in this run.":
			repeat_line = True
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
	if moffles.startswith(" "):  # Null linkers.  Make .smi compatible
		moffles = "*" + moffles + "no_mof."
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
		print("Bad MOFFLES:", moffles)
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

def cif2moffles(cif_path, intermediate_output_path=DEFAULT_OUTPUT_PATH):
	# Assemble the MOFFLES code from all of its pieces
	linkers, cat = extract_linkers(cif_path, intermediate_output_path)
	if cat is not None:
		sn_topology = extract_topology(os.path.join(intermediate_output_path, "SingleNode/topology.cgd"))
		an_topology = extract_topology(os.path.join(intermediate_output_path, "AllNode/topology.cgd"))
		if sn_topology == an_topology:
			topology = sn_topology
		else:
			topology = sn_topology + "," + an_topology
	else:
		topology = "NA"
	mof_name = os.path.splitext(os.path.basename(cif_path))[0]
	return assemble_moffles(linkers, topology, cat, mof_name=mof_name)

if __name__ == "__main__":
	args = sys.argv[1:]
	if len(args) != 1 and len(args) != 2:
		raise SyntaxError("Usage: python extract_moffles.py path_to_cif_for_analysis.cif OutputPathIfNonstandard")
	cif_file = args[0]
	
	output_systre_and_cif_path = DEFAULT_OUTPUT_PATH
	if len(args) == 2:
		output_systre_and_cif_path = args[1]
	
	print(cif2moffles(cif_file, output_systre_and_cif_path))
