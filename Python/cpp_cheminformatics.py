"""
Load cheminformatics libraries and helper utilities

Loads Open Babel, pybel, and applicable normalization wrappers

@author: Ben Bucior
"""

# Calling an external obabel binary is more expensive than using the built-in Python
# libraries but ensures that all Open Babel calls are consistent.

import sys, os
import re
from mofid.paths import openbabel_path, bin_path
# Ensure that subprocess32 is installed if running Py2
if sys.version_info[0] < 3:
	try:
		import subprocess32 as subprocess
	except:
		raise AssertionError('You must install subprocess if running Python2')
else:
	import subprocess

def runcmd(cmd_list, timeout=None):
	if timeout is None:
		return subprocess.run(cmd_list, universal_newlines=True,
			stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	else:
		return subprocess.run(cmd_list, universal_newlines=True,
			stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=timeout)

# Set up local Open Babel data environment before importing the libraries.
# CIF and other SMILES work is handled by the bin/sbu binary, called as a subprocess
os.environ['BABEL_DATADIR'] = os.path.join(openbabel_path,'data')  # directory with native EOL
TSFM_BIN = os.path.join(bin_path,'tsfm_smiles')
OBABEL_BIN = os.path.join(openbabel_path,'build','bin','obabel')

def quote(smiles_str):
	# Prepares SMILES strings for Open Babel command line calls.
	# Formerly wrapped strings within single (non-parseable) quotes, which is necessary in Bash to
	# avoid globbing/expansion but apparently not for Python's subprocess calls.  The quotes were
	# being literally passed as part of the SMILES in Linux, resulting in parsing failures.
	#return ''' + smiles_str + '''
	return smiles_str

def in_smiles(smiles_str):
	# Adds necessary quotes and Open Babel input notation for SMILES strings
	return '-:' + quote(smiles_str)

def ob_normalize(smiles):
	# Normalizes an arbitrary SMILES string with the same format and parameters as sbu.cpp
	cpp_run = runcmd([OBABEL_BIN, in_smiles(smiles), '-xi', '-ocan'])
	cpp_output = cpp_run.stdout
	if (cpp_run.stderr != '1 molecule converted\n'):
		sys.stderr.write(cpp_run.stderr + '\n')  # Re-fowarding obabel errors
	return cpp_output.rstrip()

def openbabel_replace(mol_smiles, query, replacement):
	# Perform Open Babel transforms, deletions, and/or replacements on a SMILES molecule.
	# With help from on http://baoilleach.blogspot.com/2012/08/transforming-molecules-intowellother.html
	# See also the [Daylight manual on SMARTS](http://www.daylight.com/dayhtml/doc/theory/theory.smarts.html)
	# and phmodel.cpp:208, which clarifies the possibilities of Open Babel replacements
	cpp_run = runcmd([TSFM_BIN, quote(mol_smiles), quote(query), quote(replacement)])
	cpp_output = cpp_run.stdout
	sys.stderr.write(cpp_run.stderr)  # Re-fowarding C++ errors
	return ob_normalize(cpp_output.rstrip())

def openbabel_contains(mol_smiles, query):
	# Checks if a molecule (including multi-fragment contains a SMARTS match
	cpp_run = runcmd([OBABEL_BIN, in_smiles(mol_smiles), '-s', quote(query),
		'-xi', '-ocan'])
	if (cpp_run.stderr == '1 molecule converted\n'):
		return True
	elif (cpp_run.stderr == '0 molecules converted\n'):
		return False
	else:
		sys.stderr.write(cpp_run.stderr + '\n')  # Re-fowarding obabel errors
		return False

def openbabel_formula(mol_smiles):
	# Extracts a molecular formula without relying on the pybel module
	# The .txt format prints the title: https://openbabel.org/docs/dev/FileFormats/Title_format.html
	# Note: it looks like the various --append options are in descriptors/filters.cpp, etc.
	cpp_run = runcmd([OBABEL_BIN, in_smiles(mol_smiles), '--append',
		'FORMULA', '-otxt'])
	cpp_output = cpp_run.stdout
	if (cpp_run.stderr != '1 molecule converted\n'):
		sys.stderr.write(cpp_run.stderr + '\n')  # Re-fowarding obabel errors
	return cpp_output.rstrip()

def openbabel_GetSpacedFormula(mol_smiles, delim=' '):
	# Re-implements part of OpenBabel's GetSpacedFormula method
	consolidated_formula = openbabel_formula(mol_smiles)
	split_formula = re.findall(r'[A-Z][a-z]?|d*', consolidated_formula)
	if split_formula[-1] == '':
		split_formula.pop()
	return delim.join(split_formula)
