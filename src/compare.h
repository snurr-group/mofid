#ifndef COMPARE_H
#define COMPARE_H

#include <string>

#include <openbabel/mol.h>

using namespace OpenBabel;

OBMol getMolFromCIF(const std::string& pathToFile);
bool areMolsSame(OBMol* mol1, OBMol* mol2);
bool areBondsSame(OBMol* mol1, OBMol* mol2);
bool isAtomSame(OBAtom* atom1, OBAtom* atom2);
bool isClose(double A, double B);

#endif
