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
#include <openbabel/bond.h>
#include <openbabel/math/vector3.h>

#include "config_sbu.h"
#include "compare.h"

using namespace OpenBabel;

int main(int argc, char* argv[]) 
{
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
    OBMol mol{};
    obconversion.Read(&mol);
    return mol;
}

bool areMolsSame(OBMol* mol1, OBMol* mol2) {
    bool isIdentical{true};
    if (!(mol1->NumAtoms() == mol2->NumAtoms())) {
        std::cout << "Molecular atom count mismatch" << std::endl;
        isIdentical = false;
    }
    if (mol1->GetFormula() != mol2->GetFormula()) {
        std::cout << "Molecular stoichiometric formula mismatch" << std::endl;
        isIdentical = false;
    }
    if (mol1->GetEnergy() != mol2->GetEnergy()) {
        std::cout << "Molecular energy mismatch" << std::endl;
        isIdentical = false;
    }
    if (mol1->GetMolWt() != mol2->GetMolWt()) {
        std::cout << "Molecular weight mismatch" << std::endl;
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
        return areBondsSame(mol1, mol2);
    }
}

bool areBondsSame(OBMol* mol1, OBMol* mol2) {
    bool isIdentical{true};
    if (!(mol1->NumBonds() == mol2->NumBonds())) {
        isIdentical = false;
    }
    OBBondIterator mol1Bonds{mol1->BeginBonds()};
    const OBBondIterator mol1BondsEnd{mol1->EndBonds()};
    std::map<std::set<std::vector<double>>, OBBond*> coordsToBond{};
    while (mol1Bonds != mol1BondsEnd) {
        OBBond* mol1Bond{*mol1Bonds};
        OBAtom* mol1BeginAtom{mol1Bond->GetBeginAtom()};
        OBAtom* mol1EndAtom{mol1Bond->GetEndAtom()};
        std::set<std::vector<double>> coords{{mol1BeginAtom->GetX(), mol1BeginAtom->GetY(), mol1BeginAtom->GetZ()}, {mol1EndAtom->GetX(), mol1EndAtom->GetY(), mol1EndAtom->GetZ()}};
        coordsToBond.insert(std::make_pair(coords, mol1Bond));
        ++mol1Bonds;
    }
    OBBondIterator mol2Bonds{mol2->BeginBonds()};
    const OBBondIterator mol2BondsEnd{mol2->EndBonds()};
    while (mol2Bonds != mol2BondsEnd) {
        OBBond* mol2Bond{*mol2Bonds};
        OBAtom* mol2BeginAtom{mol2Bond->GetBeginAtom()};
        OBAtom* mol2EndAtom{mol2Bond->GetEndAtom()};
        std::set<std::vector<double>> coords{{mol2BeginAtom->GetX(), mol2BeginAtom->GetY(), mol2BeginAtom->GetZ()}, {mol2EndAtom->GetX(), mol2EndAtom->GetY(), mol2EndAtom->GetZ()}};
        if (coordsToBond.count(coords) > 0) {
            OBBond* mol1Bond{coordsToBond[coords]};
            if (!isClose(mol1Bond->GetBondOrder(), mol2Bond->GetBondOrder())) {
                std::cout << "Bond order mismatch" << std::endl;
                isIdentical = false;
            }
            if (!isClose(mol1Bond->GetLength(), mol2Bond->GetLength())) {
                std::cout << "Bond length mismatch" << std::endl;
                isIdentical = false;
            }
            // std::cout << "Bond: " << "(" << mol1Bond->GetBeginAtom()->GetVector() << " " << mol1Bond->GetEndAtom()->GetVector() << ") (" << mol2Bond->GetBeginAtom()->GetVector() << " " << mol2Bond->GetEndAtom()->GetVector() << ")" << std::endl;
            if (!(isAtomSame(mol1Bond->GetBeginAtom(), mol2Bond->GetBeginAtom()) || isAtomSame(mol1Bond->GetBeginAtom(), mol2Bond->GetEndAtom()))) {
                std::cout << "Bond begin atom mismatch" << std::endl;
                isIdentical = false;
            }
            if (!(isAtomSame(mol1Bond->GetEndAtom(), mol2Bond->GetBeginAtom()) || isAtomSame(mol1Bond->GetEndAtom(), mol2Bond->GetEndAtom()))) {
                std::cout << "Bond end atom mismatch" << std::endl;
                isIdentical = false;
            }
        } else {
            std::cout << "Bond mismatch" << std::endl;
            isIdentical = false;
        }
        ++mol2Bonds;
    }
    return isIdentical;
}

bool isAtomSame(OBAtom* atom1, OBAtom* atom2) {
    if (atom1->GetAtomicNum() != atom2->GetAtomicNum()) {
        return false;
    }
    if (atom1->GetHyb() != atom2->GetHyb()) {
        return false;
    }
    if (atom1->GetIsotope() != atom2->GetIsotope()) {
        return false;
    }
    if (atom1->GetTotalDegree() != atom2->GetTotalDegree()) {
        return false;
    }
    return true;
}

bool isClose(double A, double B) {
    return (std::fabs(A - B) <= 0.0001);
}
