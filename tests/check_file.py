from math import isclose, isnan
from openbabel import pybel
import sys

if __name__ == "__main__":
    if (len(sys.argv) != 3):
        sys.exit("Usage: python check_file.py path_to_file_1 path_to_file_2")
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    if ((not file1.lower().endswith("cif")) or (not file2.lower().endswith("cif"))):
        sys.exit("Unexpected file type, expected CIF files")
    mol1 = next(pybel.readfile("cif", file1))
    mol2 = next(pybel.readfile("cif", file2))
    mismatch = False
    desc1 = mol1.calcdesc()
    desc2 = mol2.calcdesc()
    for k in desc1.keys():
        if (isnan(desc1[k]) and isnan(desc2[k])):
            continue
        if ((type(desc1[k]) is float) and (type(desc2[k]) is float)):
            if (not isclose(desc1[k], desc2[k])):
                print(f"Molecule {k} mismatch: {desc1[k]} vs. {desc2[k]}")
                mismatch = True
        elif (desc1[k] != desc2[k]):
            print(f"Molecule {k} mismatch: {desc1[k]} vs. {desc2[k]}")
            mismatch = True
    if (len(mol1.atoms) != len(mol2.atoms)):
        print(f"Atom count mismatch: {len(mol1.atoms)} vs. {len(mol2.atoms)}")
        mismatch = True
    if (mol1.calcfp().bits != mol2.calcfp().bits):
        print(f"Molecular fingerprint (FP2) mismatch: {mol1.calcfp().bits} vs. {mol2.calcfp().bits}")
        mismatch = True
    if (mol1.energy != mol2.energy):
        print(f"Molecular energy mismatch: {mol1.energy} vs. {mol2.energy}")
        mismatch = True
    if (mol1.formula != mol2.formula):
        print(f"Molecular formula mismatch: {mol1.formula} vs. {mol2.formula}")
        mismatch = True
    if (mol1.exactmass != mol2.exactmass):
        print(f"Molecular mass mismatch: {mol1.exactmass} vs. {mol2.exactmass}")
    if (mol1.molwt != mol2.molwt):
        print(f"Molecular weight mismatch: {mol1.molwt} vs. {mol2.molwt}")
    if (mol1.spin != mol2.spin):
        print(f"Spin multiplicity mismatch: {mol1.spin} vs. {mol2.spin}")
    if (mol1.charge != mol2.charge):
        print(f"Molecule charge mismatch: {mol1.charge} vs. {mol2.charge}")
    if (mismatch):
        sys.exit(1)
    else:
        sys.exit()
