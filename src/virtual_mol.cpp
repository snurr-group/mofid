// TODO: consider making the "false" returns raise an OBError/Warning
#include "virtual_mol.h"
#include "obdetails.h"
#include "framework.h"

#include <map>
#include <set>

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>

namespace OpenBabel
{

VirtualMol::VirtualMol(OBMol *parent) {
	_parent_mol = parent;
}

OBMol* VirtualMol::GetParent() {
	return _parent_mol;
}

std::set<OBAtom*> VirtualMol::GetAtoms() {
	return _atoms;
}

bool VirtualMol::HasAtom(OBAtom *a) {
	return (_atoms.find(a) != _atoms.end());
}

bool VirtualMol::AddAtom(OBAtom *a) {
	if (a->GetParent() != _parent_mol) { return false; }
	if (HasAtom(a)) { return false; }
	_atoms.insert(a);
	return true;
}

bool VirtualMol::RemoveAtom(OBAtom *a) {
	if (a->GetParent() != _parent_mol) { return false; }
	if (!HasAtom(a)) { return false; }
	_atoms.erase(a);
	return true;
}

bool VirtualMol::AddVirtualMol(VirtualMol addition) {
	if (addition.GetParent() != _parent_mol) { return false; }
	std::set<OBAtom*> atoms_to_add = addition.GetAtoms();
	for (std::set<OBAtom*>::iterator it=atoms_to_add.begin(); it!=atoms_to_add.end(); ++it) {
		_atoms.insert(*it);
	}
	return true;
}

OBMol VirtualMol::ToOBMol(bool export_bonds, bool copy_bonds) {
	// Adapts VirtualMol to an OpenBabel::OBMol
	// FIXME: does not copy paddlewheel properties.
	// Could redo formAtom/(getAtom?) with a minimalAtom if necessary.
	// Warning: the copy_bonds=false path is untested

	OBMol mol;
	typedef std::map<OBAtom*, OBAtom*> amap_t;
	amap_t virtual_to_real;
	// Copy atoms
	for (std::set<OBAtom*>::iterator it=_atoms.begin(); it!=_atoms.end(); ++it) {
		OBAtom* virtual_atom = (*it);
		OBAtom* real_atom;
		real_atom = formAtom(&mol, virtual_atom->GetVector(), virtual_atom->GetAtomicNum());
		virtual_to_real[virtual_atom] = real_atom;
	}

	if (!export_bonds) {
		return mol;  // only returning atoms
	}

	if (copy_bonds) {
		std::set<OBBond*> virtual_bonds;
		for(amap_t::iterator it=virtual_to_real.begin(); it!=virtual_to_real.end(); ++it) {
			OBAtom* v_atom = it->first;
			FOR_NBORS_OF_ATOM(n, *(v_atom)) {
				virtual_bonds.insert(v_atom->GetBond(&*n));
			}
		}
		for (std::set<OBBond*>::iterator it=virtual_bonds.begin(); it!=virtual_bonds.end(); ++it) {
			OBAtom* v_a1 = (*it)->GetBeginAtom();
			OBAtom* v_a2 = (*it)->GetEndAtom();
			int v_order = (*it)->GetBondOrder();
			formBond(&mol, virtual_to_real[v_a1], virtual_to_real[v_a2], v_order);
		}
	} else {  // recalculating bonds based on distance, etc.
		resetBonds(&mol);
	}
	return mol;
}

} // end namespace OpenBabel
