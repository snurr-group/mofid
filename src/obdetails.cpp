#include "obdetails.h"
#include "invector.h"

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/elements.h>


namespace OpenBabel
{

bool isMetal(const OBAtom* atom) {
	// Use the InChI definition of metals from https://jcheminf.springeropen.com/articles/10.1186/s13321-015-0068-4
	// Nonmetals are H, He, B, C, N, O, F, Ne, Si, P, S, Cl, Ar, Ge, As, Se, Br, Kr, Te, I, Xe, At, Rn
	const int NUM_NONMETALS = 23;
	unsigned int nonmetals[NUM_NONMETALS] = {1, 2, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18, 32, 33, 34, 35, 36, 52, 53, 54, 85, 86};

	unsigned int element = atom->GetAtomicNum();
	for (int i=0; i<NUM_NONMETALS; ++i) {
		if (nonmetals[i] == element) {
			return false;  // Atom is a nonmetal
		}
	}
	return true;  // Atom not classified as a nonmetal, so it's a metal
}

OBBond* formBond(OBMol *mol, OBAtom *begin, OBAtom *end, int order) {
	// Makes a bond between two atoms, complete with the proper accounting
	// TODO: decide how to handle cases where the bond already exists.
	// Overwrite existing bonds?  Make it, as long as it's a different periodic direction??
	OBBond* pseudo_link = NULL;
	if (!begin || !end) {  // either atom is NULL
		obErrorLog.ThrowError(__FUNCTION__, "begin or end OBAtom is undefined", obError);
	} else if (mol->GetBond(begin, end)) {
		obErrorLog.ThrowError(__FUNCTION__, "Did not generate multiply-defined bond between two atoms.", obWarning);
	} else {
		pseudo_link = mol->NewBond();
		pseudo_link->SetBegin(begin);
		pseudo_link->SetEnd(end);
		pseudo_link->SetBondOrder(order);
		// Per OBBuilder::Connect in builder.cpp, we need to also update the atoms' bonding.
		// Otherwise, our OBAtomAtomIter will not operate properly (bonds will not propagate back to the atoms).
		begin->AddBond(pseudo_link);
		end->AddBond(pseudo_link);
	}
	return pseudo_link;
}

OBAtom* formAtom(OBMol *mol, vector3 loc, int element) {
	// Makes a new atom with a specified location and atomic number
	OBAtom* atom = mol->NewAtom();
	atom->SetVector(loc);
	atom->SetAtomicNum(element);
	atom->SetType(OBElements::GetName(element));
	return atom;
}

int deleteBonds(OBMol *mol, bool only_metals) {
	// Deletes bonds in a molecule, optionally only deleting bonds to metal atoms.
	// Returns the number of bond deletions.
	std::vector<OBMol> fragments;
	obErrorLog.ThrowError(__FUNCTION__, "Bond deletion enabled", obDebug);
	// Same spirit as OBMol::DeleteHydrogens()
	std::vector<OBBond*> delbonds;
	FOR_ATOMS_OF_MOL(a, *mol) {
		if (isMetal(&*a) || !only_metals) {
			FOR_BONDS_OF_ATOM(b, *a) {
				if (!inVector(&*b, delbonds)) {
					// avoid double deletion in the case of metal-metal bonds
					delbonds.push_back(&*b);
				}
			}
		}
	}
	mol->BeginModify();
	std::stringstream deletionMsg;
	for (std::vector<OBBond*>::iterator itb = delbonds.begin(); itb != delbonds.end(); ++itb) {
		deletionMsg << "Deleting bond with atomic numbers ";
		deletionMsg << (*itb)->GetBeginAtom()->GetAtomicNum() << " and ";
		deletionMsg << (*itb)->GetEndAtom()->GetAtomicNum();
		deletionMsg << " (OBBond * " << static_cast<void *>(&*itb) << ")\n";
		mol->DeleteBond(*itb);  // We should also consider preserving the connection point instead of outright deleting the bond
	}
	obErrorLog.ThrowError(__FUNCTION__, deletionMsg.str(), obDebug);
	mol->EndModify();
	return delbonds.size();
}

} // end namespace OpenBabel
