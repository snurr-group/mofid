#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Digest MOFs into their nodes and linkers

Use the InChI/CSD bonding criteria for now, ignoring bonds to metals.
Later, Use the TOPOS bonding criteria for inter/intramolecular bonding to
decompose MOF structures into their building blocks.  Eventually, this code
should also find the MOF topology as a user-friendly classification tool.

@author: Ben Bucior
"""

import pybel
import openbabel

def get_unique_smiles(filename):
    # Return a list of unique canonical SMILES representations of a molecule
    mol_class = pybel.readfile("cif", filename).next()
    mol_class = mol_class.OBMol
    separate_mols = mol_class.Separate()
    unique_smiles = []
    for molecule in separate_mols:
        #print molecule.AddHydrogens()
        #for atom in openbabel.OBMolAtomIter(molecule):
            #print atom.GetAtomicNum(), atom.GetSpinMultiplicity()
            #atom.ForceImplH()
        mol_smiles = pybel.Molecule(molecule).write("can").split("\t")[0]
        if mol_smiles not in unique_smiles:
            unique_smiles.append(mol_smiles)
    return unique_smiles


import os
print os.getcwd()
CIF_PATH = "../../Data/RingCIFs/"
print "HKUST-1", get_unique_smiles(os.path.join(CIF_PATH, "DOTSOV_clean.cif"))
print "IRMOF-1", get_unique_smiles(os.path.join(CIF_PATH, "P1-IRMOF-1.cif"))
print "IRMOF-1", get_unique_smiles(os.path.join(CIF_PATH, "VUSKEA_clean.cif"))
print "Pillar (noncatenated EDOMAM)", get_unique_smiles(os.path.join(CIF_PATH, "hypotheticalMOF_5071180_i_2_j_25_k_1_m_0.cif"))
