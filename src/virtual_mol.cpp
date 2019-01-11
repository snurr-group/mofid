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

OBMol VirtualMol::ToOBMol(bool export_bonds, bool copy_bonds) {
	// Copies VirtualMol to an OpenBabel::OBMol
	// WARNING: the copy_bonds=false path is untested

	OBMol mol = initMOFwithUC(_parent_mol);
	typedef std::map<OBAtom*, OBAtom*> amap_t;
	amap_t virtual_to_mol;
	// Copy atoms
	for (std::set<OBAtom*>::iterator it=_atoms.begin(); it!=_atoms.end(); ++it) {
		OBAtom* virtual_atom = (*it);
		OBAtom* mol_atom;
		mol_atom = formAtom(&mol, virtual_atom->GetVector(), virtual_atom->GetAtomicNum());
		virtual_to_mol[virtual_atom] = mol_atom;

		// Also copy paddlewheel detection.  Based on framework.cpp:resetBonds
		bool is_paddlewheel = virtual_atom->HasData("Paddlewheel");
		if (is_paddlewheel) {
			OBPairData *dp = new OBPairData;
			dp->SetAttribute("Paddlewheel");
			mol_atom->SetData(dp);
		}
	}

	if (!export_bonds) {
		return mol;  // only returning atoms
	}

	if (copy_bonds) {
		std::set<OBBond*> virtual_bonds;
		for(amap_t::iterator it=virtual_to_mol.begin(); it!=virtual_to_mol.end(); ++it) {
			OBAtom* v_atom = it->first;
			// Only export bonds if both atoms (v_atom and n) are in the VirtualMol
			FOR_NBORS_OF_ATOM(n, *(v_atom)) {
				if (HasAtom(&*n)) {
					virtual_bonds.insert(v_atom->GetBond(&*n));
				}
			}
		}
		for (std::set<OBBond*>::iterator it=virtual_bonds.begin(); it!=virtual_bonds.end(); ++it) {
			OBAtom* v_a1 = (*it)->GetBeginAtom();
			OBAtom* v_a2 = (*it)->GetEndAtom();
			int v_order = (*it)->GetBondOrder();
			formBond(&mol, virtual_to_mol[v_a1], virtual_to_mol[v_a2], v_order);
		}
	} else {  // recalculating bonds based on distance, etc.
		resetBonds(&mol);
	}
	return mol;
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
