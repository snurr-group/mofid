#ifndef COMPARE_H
#define COMPARE_H

#include <vector>
#include <string>

#include <openbabel/mol.h>
#include <openbabel/bond.h>

using namespace OpenBabel;

OBMol getMolFromCIF(const std::string& pathToFile);
bool areMolsSame(OBMol* mol1, OBMol* mol2);
bool areAtomsSame(OBMol* mol1, OBMol* mol2, int precision);
bool isAtomSame(OBMol* mol1, OBAtom* atom1, OBMol* mol2, OBAtom* atom2);
std::vector<OBBond*> getBonds(OBMol* mol, OBAtom* atom);
bool areBondsSame(OBMol* mol1, OBAtom* atom1, OBMol* mol2, OBAtom* atom2);
void printBonds(OBMol* mol, const std::string& name, int precision);
void printBond(OBBond* bond, const std::string& name, int precision);
bool isClose(double A, double B);

#endif
