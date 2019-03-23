#include <vector>
#include <map>

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
	changeAtomElement(atom, element);
	return atom;
}

void changeAtomElement(OBAtom* atom, int element) {
	// Changes an atom's atomic number
	if (atom == NULL) { obErrorLog.ThrowError(__FUNCTION__, "Skipping invalid atom", obWarning); }
	atom->SetAtomicNum(element);
	atom->SetType(OBElements::GetName(element));
}

int deleteBonds(OBMol *mol, bool only_metals) {
	// TODO: consider refactoring another method as std::vector< std::pair<OBAtom*,OBAtom*> >,
	// and use deletebonds(...).size() to get the int version
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

bool subtractMols(OBMol *mol, OBMol *subtracted) {
	// Subtracts all atoms from a second molecule present in the first
	// Returns true if successful, false if failure (extra atoms, etc)
	// Modifies *mol if true, reverts back to original *mol if false
	// This code assumes that the molecule is well-defined (no multiply-defined atoms)

	std::vector<OBAtom*> deleted;  // Consider implementing as a queue
	FOR_ATOMS_OF_MOL(a2, *subtracted) {
		OBAtom* match = atomInOtherMol(&*a2, mol);
		if (!match) {  // *subtracted has an extra atom
			obErrorLog.ThrowError(__FUNCTION__, "Attempted to subtract an atom not present in the original molecule", obWarning);
			return false;  // *mol not modified
		}
		deleted.push_back(&*match);
	}

	// It's okay to delete the atoms, since they should hold different (cloned) pointers
	mol->BeginModify();
	for (std::vector<OBAtom*>::iterator it = deleted.begin(); it != deleted.end(); ++it) {
		mol->DeleteAtom(*it);
	}
	mol->EndModify();
	return true;
}

bool atomsEqual(const OBAtom &atom1, const OBAtom &atom2) {
	// Loosely based on the data copied during OBAtom::Duplicate
	// Atoms are the same if these properties (excluding bond-related) are equivalent
	return atom1.GetAtomicNum() == atom2.GetAtomicNum() &&
			atom1.GetVector().IsApprox(atom2.GetVector(), 1.0e-6) &&
			atom1.GetIsotope() == atom2.GetIsotope();
}

OBAtom* atomInOtherMol(OBAtom *atom, OBMol *mol) {
	// Is an atom within a molecule?  (as defined by atomsEqual conditions)
	// If so, return the pointer of the copycat.  Otherwise, return NULL
	// This has the nice benefit that the usage is also compatible with bool
	// ( if (atomInOtherMol(atomp, molp)) {}  will work, too.)
	FOR_ATOMS_OF_MOL(a, *mol) {
		if (atomsEqual(*a, *atom)) {
			return &*a;
		}
	}
	return NULL;
}

bool isSubMol(OBMol *sub, OBMol *super) {
	// Are the atoms in *sub a subset of the atoms in *super?
	FOR_ATOMS_OF_MOL(a, *sub) {
		if (!atomInOtherMol(&*a, super)) {  // NULL if atom not found
			return false;
		}
	}
	return true;
}

std::map<int,int> getNumericFormula(OBMol *mol) {
	// What is the molecular formula, based on atomic numbers?
	// Returns <atomic number of an element, atom count in *mol>
	std::map<int,int> formula;
	FOR_ATOMS_OF_MOL(a, *mol) {
		int element = a->GetAtomicNum();
		if (formula.find(element) == formula.end()) {
			formula[element] = 0;
		}
		formula[element] += 1;
	}
	return formula;
}

std::string rtrimWhiteSpace(const std::string str) {
	// Right-trims white space from a string, per obconversion.cpp and consensus from SO
	std::string trimmed(str);
	std::string::size_type notwhite = trimmed.find_last_not_of(" \t\n\r");
	trimmed.erase(notwhite+1);
	return trimmed;
}

} // end namespace OpenBabel
