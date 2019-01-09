#include "pseudo_atom.h"
#include "virtual_mol.h"

#include <map>

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>

namespace OpenBabel
{

PseudoAtomMap::PseudoAtomMap(OBMol *pseudo, OBMol *orig) {
	_pseudo_mol = pseudo;
	_full_mol = orig;
}

VirtualMol& PseudoAtomMap::operator[] (PseudoAtom i) {
	return _mapping[i];
}

OBMol PseudoAtomMap::ToCombinedMol(bool export_bonds, bool copy_bonds) {
	VirtualMol combined(_full_mol);
	for (std::map<PseudoAtom, VirtualMol>::iterator it=_mapping.begin(); it!=_mapping.end(); ++it) {
		combined.AddVirtualMol(it->second);
	}
	return combined.ToOBMol(export_bonds, copy_bonds);
}

void PseudoAtomMap::RemoveAtom(PseudoAtom atom) {
	_mapping.erase(atom);
}

} // end namespace OpenBabel
