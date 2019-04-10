"""
Load cheminformatics libraries and helper utilities

Loads Open Babel, pybel, and applicable normalization wrappers

@author: Ben Bucior
"""

import os
from mofid.paths import openbabel_path
# Set up local Open Babel data environment before importing the libraries.
# CIF and other SMILES work is handled by the bin/sbu binary, called as a subprocess
os.environ['BABEL_DATADIR'] = os.path.join(openbabel_path,'data')  # directory with native EOL
import pybel  # Read SMILES to calculate molecular formulas, run SMARTS-based OBChemTsfm, etc.

def ob_normalize(smiles):
	# Normalizes an arbitrary SMILES string with the same format and parameters as sbu.cpp
	ob_mol = pybel.readstring('smi', smiles)
	return ob_mol.write('can', opt={'i': True}).rstrip()

def openbabel_replace(mol_smiles, query, replacement):
	# Perform Open Babel transforms, deletions, and/or replacements on a SMILES molecule.
	# With help from on http://baoilleach.blogspot.com/2012/08/transforming-molecules-intowellother.html
	# See also the [Daylight manual on SMARTS](http://www.daylight.com/dayhtml/doc/theory/theory.smarts.html)
	# and phmodel.cpp:208, which clarifies the possibilities of Open Babel replacements
	transform = pybel.ob.OBChemTsfm()
	success = transform.Init(query, replacement)
	assert success
	mol = pybel.readstring('smi', mol_smiles)
	transform.Apply(mol.OBMol)
	mol.OBMol.UnsetAromaticPerceived()
	mol.OBMol.Kekulize()
	return ob_normalize(mol.write('can'))

def openbabel_contains(mol_smiles, query):
	# Checks if a molecule (including multi-fragment contains a SMARTS match
	matcher = pybel.ob.OBSmartsPattern()
	success = matcher.Init(query)
	assert success
	mol = pybel.readstring('smi', mol_smiles)
	return matcher.Match(mol.OBMol)

