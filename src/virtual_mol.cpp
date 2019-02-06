#include "virtual_mol.h"
#include "obdetails.h"
#include "framework.h"

#include <vector>
#include <map>
#include <set>
#include <utility>  // std::pair
#include <stack>

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

VirtualMol::VirtualMol(OBAtom *single_atom) {
	_parent_mol = single_atom->GetParent();
	AddAtom(single_atom);
}

OBMol* VirtualMol::GetParent() {
	return _parent_mol;
}

int VirtualMol::NumAtoms() {
	return _atoms.size();
}

std::set<OBAtom*> VirtualMol::GetAtoms() {
	return _atoms;
}

bool VirtualMol::HasAtom(OBAtom *a) {
	return (_atoms.find(a) != _atoms.end());
}

bool VirtualMol::AddAtom(OBAtom *a) {
	if (a->GetParent() != _parent_mol) {
		obErrorLog.ThrowError(__FUNCTION__, "Atom does not have the same parent OBMol as the VirtualMol.", obError);
		return false;
	}
	if (HasAtom(a)) { return false; }
	_atoms.insert(a);
	return true;
}

bool VirtualMol::RemoveAtom(OBAtom *a) {
	if (a->GetParent() != _parent_mol) {
		obErrorLog.ThrowError(__FUNCTION__, "Atom does not have the same parent OBMol as the VirtualMol.", obError);
		return false;
	}
	if (!HasAtom(a)) {
		obErrorLog.ThrowError(__FUNCTION__, "Tried to delete a nonexistent atom from VirtualMol.", obError);
		return false;
	}
	_atoms.erase(a);
	return true;
}

bool VirtualMol::AddVirtualMol(VirtualMol addition) {
	if (addition.GetParent() != _parent_mol) {
		obErrorLog.ThrowError(__FUNCTION__, "VirtualMol parents do not match", obWarning);
		return false;
	}
	std::set<OBAtom*> atoms_to_add = addition.GetAtoms();
	for (std::set<OBAtom*>::iterator it=atoms_to_add.begin(); it!=atoms_to_add.end(); ++it) {
		_atoms.insert(*it);
	}
	return true;
}

int VirtualMol::ImportCopiedFragment(OBMol *fragment) {
	// Adds an OBMol based on the former methodology of searching for atoms with the same
	// element, position, etc.
	std::vector<OBAtom*> atoms_to_add;
	atoms_to_add.reserve(fragment->NumAtoms());
	FOR_ATOMS_OF_MOL(a, *fragment) {
		OBAtom* parent_a;
		parent_a = atomInOtherMol(&*a, _parent_mol);
		if (!parent_a) {  // NULL: atom does not exist
			obErrorLog.ThrowError(__FUNCTION__, "Cannot import fragment with missing atom.", obError);
			return 0;  // error, and no modifications to VirtualMol
		}
		atoms_to_add.push_back(parent_a);
	}
	for (std::vector<OBAtom*>::iterator it=atoms_to_add.begin(); it!=atoms_to_add.end(); ++it) {
		AddAtom(*it);
	}
	return atoms_to_add.size();
}

ConnIntToExt VirtualMol::GetExternalBondsOrConns() {
	// Provides the mapping of VirtualMol atoms bonded to nearest-neighbor external atoms in the parent molecule.
	// This function replaces getLinksToExt from sbu.cpp.
	// WARNING: if this function is run on a simplified_net, it will consider connection sites as
	// external unless they're part of the VirtualMol
	ConnIntToExt connections;
	for (std::set<OBAtom*>::iterator it=_atoms.begin(); it!=_atoms.end(); ++it) {
		FOR_NBORS_OF_ATOM(nbor, *it) {
			if (!HasAtom(&*nbor)) {
				std::pair<OBAtom*, OBAtom*> bond(*it, &*nbor);
				connections.insert(bond);
			}
		}
	}
	return connections;
}

void VirtualMol::CopyToMappedMol(MappedMol *dest, bool export_bonds, bool copy_bonds) {
	// Copies VirtualMol to a destination MappedMol
	// WARNING: the copy_bonds=false path is untested

	if (dest->origin_molp) {
		obErrorLog.ThrowError(__FUNCTION__, "Overwriting contents of existing MappedMol", obWarning);
	}
	dest->origin_molp = _parent_mol;
	dest->mol_copy = initMOFwithUC(_parent_mol);
	OBMol* pmol_copied = &(dest->mol_copy);  // for convenience
	dest->origin_to_copy.clear();
	dest->copy_to_origin.clear();
	dest->copy_pa_to_multiple.clear();

	// Copy atoms
	for (std::set<OBAtom*>::iterator it=_atoms.begin(); it!=_atoms.end(); ++it) {
		OBAtom* virtual_atom = (*it);
		OBAtom* mol_atom = formAtom(pmol_copied, virtual_atom->GetVector(), virtual_atom->GetAtomicNum());
		dest->origin_to_copy[virtual_atom] = mol_atom;
		dest->copy_to_origin[mol_atom] = virtual_atom;
		dest->copy_pa_to_multiple[mol_atom] = VirtualMol(virtual_atom);

		// Also copy paddlewheel detection.  Based on framework.cpp:resetBonds
		bool is_paddlewheel = virtual_atom->HasData("Paddlewheel");
		if (is_paddlewheel) {
			OBPairData *dp = new OBPairData;
			dp->SetAttribute("Paddlewheel");
			mol_atom->SetData(dp);
		}
	}

	if (!export_bonds) {
		return;  // Done: we've already copied the atoms
	}

	if (copy_bonds) {
		FOR_BONDS_OF_MOL(b, *_parent_mol) {
			OBAtom* v_a1 = b->GetBeginAtom();
			OBAtom* v_a2 = b->GetEndAtom();
			int v_order = b->GetBondOrder();
			if (HasAtom(v_a1) && HasAtom(v_a2)) {
				OBAtom* copied_a1 = dest->origin_to_copy[v_a1];
				OBAtom* copied_a2 = dest->origin_to_copy[v_a2];
				formBond(pmol_copied, copied_a1, copied_a2, v_order);
			}
		}
	} else {  // recalculating bonds based on distance, etc.
		resetBonds(pmol_copied);
	}
	return;
}

OBMol VirtualMol::ToOBMol(bool export_bonds, bool copy_bonds) {
	// Wrapper for CopyToMappedMol() if we're only interested in an unmapped OBMol copy
	MappedMol temp_map;
	CopyToMappedMol(&temp_map, export_bonds, copy_bonds);
	return temp_map.mol_copy;
}

void VirtualMol::ToCIF(const std::string &filename, bool write_bonds) {
	// Writes a CIF file from ToOBMol at the specified file path
	OBMol mol_for_export = ToOBMol(write_bonds);
	writeCIF(&mol_for_export, filename, write_bonds);
}

std::vector<VirtualMol> VirtualMol::Separate() {
	// Functions like OBMol::Separate, giving a vector of distinct, unconnected molecular fragments
	std::vector<VirtualMol> fragments;  // return value

	// Initalize visited: none visited at the beginning
	std::map<OBAtom*, bool> visited;
	for (AtomSet::iterator it=_atoms.begin(); it!=_atoms.end(); ++it) {
		visited[*it] = false;
	}

	// Outer loop through the full list of atoms to start new fragments
	for (std::map<OBAtom*, bool>::iterator next_atom=visited.begin(); next_atom!=visited.end(); ++next_atom) {
		if (next_atom->second) { continue; }  // not a new fragment: already visited
		VirtualMol curr_fragment(GetParent());

		// DFS through the network to map out the fragment
		std::stack<OBAtom*> to_visit;
		to_visit.push(next_atom->first);
		visited[next_atom->first] = true;

		while (!to_visit.empty()) {
			OBAtom* curr_atom = to_visit.top();
			to_visit.pop();
			curr_fragment.AddAtom(curr_atom);

			FOR_NBORS_OF_ATOM(nbor, *curr_atom) {
				// Make sure we're not iterating on external parent_mol atoms
				// or creating an infinite loop by continuously visiting the same atoms.
				if (this->HasAtom(&*nbor) && !visited[&*nbor]) {
					visited[&*nbor] = true;
					to_visit.push(&*nbor);
				}
			}
		}
		fragments.push_back(curr_fragment);
	}

	return fragments;
}

} // end namespace OpenBabel
