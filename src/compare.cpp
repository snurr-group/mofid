#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <utility>
#include <map>
#include <set>

#include <openbabel/obconversion.h>
#include <openbabel/generic.h>
#include <openbabel/atom.h>
#include <openbabel/math/vector3.h>

#include "config_sbu.h"
#include "compare.h"

using namespace OpenBabel;

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Incorrect number of arguments. Need to specify two CIFs." << std::endl;
        return 2;
    }
#ifdef _WIN32
	_putenv_s("BABEL_DATADIR", LOCAL_OB_DATADIR);
	_putenv_s("BABEL_LIBDIR", LOCAL_OB_LIBDIR);
#else
	setenv("BABEL_DATADIR", LOCAL_OB_DATADIR, 1);
	setenv("BABEL_LIBDIR", LOCAL_OB_LIBDIR, 1);
#endif
    OBMol mol1{getMolFromCIF(argv[1])};
    OBMol mol2{getMolFromCIF(argv[2])};
    if (areMolsSame(&mol1, &mol2)) {
        return 0;
    } else {
        return 1;
    }
}

OBMol getMolFromCIF(const std::string& pathToFile) {
    std::fstream cif{pathToFile};
    OBConversion obconversion(&cif);
    obconversion.SetInFormat("CIF");
    obconversion.SetOptions("Bbs", OBConversion::INOPTIONS);
    OBMol mol{};
    obconversion.Read(&mol);
    return mol;
}

bool areMolsSame(OBMol* mol1, OBMol* mol2) {
    bool isIdentical{true};
    if (!(mol1->NumAtoms() == mol2->NumAtoms())) {
        std::cout << "Molecular atom count mismatch" << " mol1: " << mol1->NumAtoms() << " mol2: " << mol2->NumAtoms() << std::endl;
        isIdentical = false;
    }
    if (!(mol1->NumBonds() == mol2->NumBonds())) {
        std::cout << "Molecular bond count mismatch" << " mol1: " << mol1->NumBonds() << " mol2: " << mol2->NumBonds() << std::endl;
        isIdentical = false;
    }
    if (mol1->GetFormula() != mol2->GetFormula()) {
        std::cout << "Molecular stoichiometric formula mismatch" << " " << mol1->GetFormula() << " " << mol2->GetFormula() << std::endl;
        isIdentical = false;
    }
    OBUnitCell* mol1UC{static_cast<OBUnitCell*>(mol1->GetData(OBGenericDataType::UnitCell))};
    OBUnitCell* mol2UC{static_cast<OBUnitCell*>(mol2->GetData(OBGenericDataType::UnitCell))};
    if (mol1UC->GetLatticeType() != mol2UC->GetLatticeType()) {
        std::cout << "Molecular unit cell lattice type mismatch" << std::endl;
        isIdentical = false;
    }
    if (!isIdentical) {
        return false;
    } else{
        return areAtomsSame(mol1, mol2, 3);
    }
}

bool areAtomsSame(OBMol* mol1, OBMol* mol2, int precision) {
    OBAtomIterator mol1Atoms{mol1->BeginAtoms()};
    const OBAtomIterator mol1AtomsEnd{mol1->EndAtoms()};
    std::map<unsigned int, std::set<OBAtom*>> atomicNumToAtoms{};
    while (mol1Atoms != mol1AtomsEnd) {
        OBAtom* mol1Atom{*mol1Atoms};
        unsigned int atomicNum{mol1Atom->GetAtomicNum()};
        if (auto search = atomicNumToAtoms.find(atomicNum); search != atomicNumToAtoms.end()) {
            atomicNumToAtoms[atomicNum].insert(mol1Atom);
        } else {
            std::set<OBAtom*> atoms{mol1Atom};
            atomicNumToAtoms.insert(std::make_pair(atomicNum, atoms));
        }
        ++mol1Atoms;
    }
    OBAtomIterator mol2Atoms{mol2->BeginAtoms()};
    const OBAtomIterator mol2AtomsEnd{mol2->EndAtoms()};
    while (mol2Atoms != mol2AtomsEnd) {
        OBAtom* mol2Atom{*mol2Atoms};
        unsigned int atomicNum{mol2Atom->GetAtomicNum()};
        if (auto search = atomicNumToAtoms.find(atomicNum); search != atomicNumToAtoms.end()) {
            bool foundMatch{false};
            for (OBAtom* mol1Atom : search->second) {
                if (isAtomSame(mol1, mol1Atom, mol2, mol2Atom)) {
                    foundMatch = true;
                    search->second.erase(mol1Atom);
                    break;
                }
            }
            if (!foundMatch) {
                std::cout << "Atom match not found " << mol2Atom->GetAtomicNum() << std::endl;
            }
        } else {
            std::cout << "Atomic number " << atomicNum << " not found" << std::endl;
        }
        ++mol2Atoms;
    }
    return true;
}

bool isAtomSame(OBMol* mol1, OBAtom* atom1, OBMol* mol2, OBAtom* atom2) {
    if (atom1->GetAtomicNum() != atom2->GetAtomicNum()) {
        return false;
    }
    if (atom1->GetHyb() != atom2->GetHyb()) {
        return false;
    }
    if (atom1->GetIsotope() != atom2->GetIsotope()) {
        return false;
    }
    if (!areBondsSame(mol1, atom1, mol2, atom2)) {
        return false;
    }
    return true;
}

std::vector<OBBond*> getBonds(OBMol* mol, OBAtom* atom) {
    OBBondIterator molBonds{mol->BeginBonds()};
    const OBBondIterator molBondsEnd{mol->EndBonds()};
    std::vector<OBBond*> bonds{};
    while (molBonds != molBondsEnd) {
        OBBond* molBond{*molBonds};
        if ((molBond->GetBeginAtom() == atom) || (molBond->GetEndAtom() == atom)) {
            bonds.push_back(molBond);
        }
        ++molBonds;
    }
    return bonds;
}

bool areBondsSame(OBMol* mol1, OBAtom* atom1, OBMol* mol2, OBAtom* atom2) {
    std::map<unsigned int, int> atomicNumToCount{};
    for (OBBond* atom1Bond : getBonds(mol1, atom1)) {
        OBAtom* endAtom{};
        if (atom1Bond->GetBeginAtom() != atom1) {
            endAtom = atom1Bond->GetBeginAtom();
        } else {
            endAtom = atom1Bond->GetEndAtom();
        }
        unsigned int atomicNum{endAtom->GetAtomicNum()};
        if (auto search = atomicNumToCount.find(atomicNum); search != atomicNumToCount.end()) {
            ++atomicNumToCount[atomicNum];
        } else {
            atomicNumToCount.insert(std::make_pair(atomicNum, 1));
        }
    }
    for (OBBond* atom2Bond : getBonds(mol2, atom2)) {
        OBAtom* endAtom{};
        if (atom2Bond->GetBeginAtom() != atom2) {
            endAtom = atom2Bond->GetBeginAtom();
        } else {
            endAtom = atom2Bond->GetEndAtom();
        }
        unsigned int atomicNum{endAtom->GetAtomicNum()};
        if (auto search = atomicNumToCount.find(atomicNum); search != atomicNumToCount.end()) {
            --atomicNumToCount[atomicNum];
            if (atomicNumToCount[atomicNum] == 0) {
                atomicNumToCount.erase(atomicNum);
            }
        } else {
            return false;
        }
    }
    return true;
}

void printBonds(OBMol* mol, const std::string& name, int precision) {
    OBBondIterator molBonds{mol->BeginBonds()};
    const OBBondIterator molBondsEnd{mol->EndBonds()};
    while (molBonds != molBondsEnd) {
        printBond(*molBonds, name, precision); 
        ++molBonds;
    }
}

void printBond(OBBond* bond, const std::string& name, int precision) {
        OBAtom* beginAtom{bond->GetBeginAtom()};
        int beginAtomX{static_cast<int>(std::trunc(beginAtom->GetX() * std::pow(10, precision)))};
        int beginAtomY{static_cast<int>(std::trunc(beginAtom->GetY() * std::pow(10, precision)))};
        int beginAtomZ{static_cast<int>(std::trunc(beginAtom->GetZ() * std::pow(10, precision)))};
        OBAtom* endAtom{bond->GetEndAtom()};
        int endAtomX{static_cast<int>(std::trunc(endAtom->GetX() * std::pow(10, precision)))};
        int endAtomY{static_cast<int>(std::trunc(endAtom->GetY() * std::pow(10, precision)))};
        int endAtomZ{static_cast<int>(std::trunc(endAtom->GetZ() * std::pow(10, precision)))};
        std::cout << name << ": (" << beginAtomX << " " << beginAtomY << " " << beginAtomZ << ") (" << endAtomX << " " << endAtomY << " " << endAtomZ << ") order: " << bond->GetBondOrder() << std::endl;
}

bool isClose(double A, double B) {
    return (std::fabs(A - B) <= 0.0001);
}
