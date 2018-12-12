#include "obdetails.h"

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/elements.h>


namespace OpenBabel
{

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


}
