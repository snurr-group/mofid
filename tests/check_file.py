import logging
import sys
from math import isclose, isnan
from openbabel import openbabel
from openbabel import pybel


class CIFComparer:
    def __init__(self, path_to_file1, path_to_file2):
        logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.WARNING)
        self.obmol1 = self._get_OBMol_from_CIF(path_to_file1)
        self.obmol2 = self._get_OBMol_from_CIF(path_to_file2)
        self.mol1 = self._get_mol_from_CIF(path_to_file1)
        self.mol2 = self._get_mol_from_CIF(path_to_file2)
        self.mismatch = False

    def compare(self):
        self._compare_descriptors()
        self._compare_atoms()
        self._compare_bonds()
        self._compare_molecules()
        if (self.mismatch):
            sys.exit(1)
        else:
            sys.exit()

    def _get_OBMol_from_CIF(self, path_to_file):
        ob_conversion = openbabel.OBConversion()
        ob_conversion.SetInAndOutFormats("cif", "mol")
        mol = openbabel.OBMol()
        ob_conversion.ReadFile(mol, path_to_file)
        return mol

    def _get_mol_from_CIF(self, path_to_file):
        if ((not path_to_file1.lower().endswith("cif")) or (not path_to_file2.lower().endswith("cif"))):
            sys.exit("Unexpected file type, expected CIF files")
        return next(pybel.readfile("cif", path_to_file))

    def _compare_descriptors(self):
        desc1 = self.mol1.calcdesc()
        desc2 = self.mol2.calcdesc()
        for k in desc1.keys():
            if (isnan(desc1[k]) and isnan(desc2[k])):
                continue
            if ((type(desc1[k]) is float) and (type(desc2[k]) is float)):
                if (not isclose(desc1[k], desc2[k])):
                    logging.warning(f"molecule %s self.mismatch, %s vs. %s", k, desc1[k], desc2[k])
                    self.mismatch = True
            elif (desc1[k] != desc2[k]):
                logging.warning(f"molecule %s self.mismatch, %s vs. %s", k, desc1[k], desc2[k])
                self.mismatch = True

    def _compare_atoms(self):
        if (len(self.mol1.atoms) != len(self.mol2.atoms)):
            logging.warning("atom count self.mismatch, %s vs. %s", len(self.mol1.atoms), len(self.mol2.atoms))
            self.mismatch = True
        else:
            for i in range(len(self.mol1.atoms)):
                atom1 = self.mol1.atoms[i]
                atom2 = self.mol2.atoms[i]
                for (attr1, attr2, message) in zip((atom1.atomicnum, atom1.exactmass, atom1.formalcharge, atom1.heavydegree, atom1.heterodegree, atom1.hyb, atom1.isotope, atom1.partialcharge, atom1.spin, atom1.type, atom1.degree), (atom2.atomicnum, atom2.exactmass, atom2.formalcharge, atom2.heavydegree, atom2.heterodegree, atom2.hyb, atom2.isotope, atom2.partialcharge, atom2.spin, atom2.type, atom2.degree), ("number", "exact mass", "formal charge", "attached non-hydrogens", "attached heteroatoms", "hybridization", "isotope", "partial charge", "spin", "type", "explicit connections")):
                    if (attr1 != attr2):
                        logging.warning("atom %s %s mismatch, %s vs. %s", i, message, attr1, attr2)
                        self.mismatch = True

    def _compare_bonds(self):
        i = 0
        for (bond1, bond2) in zip(openbabel.OBMolBondIter(self.obmol1), openbabel.OBMolBondIter(self.obmol2)):
            i += 1
            for (attr1, attr2, message) in zip((bond1.GetBondOrder(), round(bond1.GetLength(), 9), bond1.GetBeginAtom().GetAtomicNum(), bond1.GetEndAtom().GetAtomicNum()), (bond2.GetBondOrder(), round(bond2.GetLength(), 9), bond2.GetBeginAtom().GetAtomicNum(), bond2.GetEndAtom().GetAtomicNum()), ("order", "length", "begin atom", "end atom")):
                if (attr1 != attr2):
                    logging.warning("bond %s %s mismatch, %s vs. %s", i, message, attr1, attr2)
                    self.mismatch = True

    def _compare_molecules(self):
        for (attr1, attr2, message) in zip((self.mol1.calcfp().bits, self.mol1.energy, self.mol1.formula, self.mol1.exactmass, self.mol1.molwt, self.mol1.spin, self.mol1.charge), (self.mol2.calcfp().bits, self.mol2.energy, self.mol2.formula, self.mol2.exactmass, self.mol2.molwt, self.mol2.spin, self.mol2.charge), ("fingerprint (FP2)", "energy", "formula", "mass", "weight", "spin multiplicity", "charge")):
            if (attr1 != attr2):
                logging.warning("atom %s %s mismatch, %s vs. %s", i, message, attr1, attr2)
                self.mismatch = True

if __name__ == "__main__":
    if (len(sys.argv) != 3):
        sys.exit("Usage: python check_file.py path_to_file_1 path_to_file_2")
    path_to_file1 = sys.argv[1]
    path_to_file2 = sys.argv[2]
    comparer = CIFComparer(path_to_file1, path_to_file2)
    comparer.compare()
