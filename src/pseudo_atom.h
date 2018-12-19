/**********************************************************************
pseudo_atom.h - Pseudo atoms in a simplified MOF topology
***********************************************************************/

#ifndef PSEUDO_ATOM_H
#define PSEUDO_ATOM_H

#include <map>

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>

#include "virtual_mol.h"

namespace OpenBabel
{
// forward declarations
class OBAtom;

typedef OBAtom* PseudoAtom;

class PseudoAtomMap {
// Maps PseudoAtoms in a simplified _pseudo_mol to sets of atoms in the original _full_mol
private:
	OBMol *_pseudo_mol;
	OBMol *_full_mol;
	std::map<PseudoAtom, VirtualMol> _mapping;  // between _pseudo_mol and _full_mol
public:
	// PseudoAtomMap() = delete;  // this is difficult to work with.  Just set to NULL by default
	PseudoAtomMap(OBMol *psuedo = NULL, OBMol *orig = NULL);
	// void MapOneToOne();  // This will be implemented when the simplified net is created
	OBMol ToCombinedMol(bool export_bonds = true, bool copy_bonds = true);
	VirtualMol& operator[] (PseudoAtom i);
	// Possibly also a helper utility to reassign atom classifications in _full_mol.
	// It would move an original OBAtom* from one <PA,VirtualMol> to a different one.
	// A utility to find unique SMILES may also be useful, but something to consider later.
};

} // end namespace OpenBabel
#endif // PSEUDO_ATOM_H

//! \file pseudo_atom.h
//! \brief pseudo_atom.h - Pseudo atoms in a simplified MOF topology
