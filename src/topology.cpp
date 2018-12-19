#include "topology.h"
#include "framework.h"
#include "virtual_mol.h"
#include "pseudo_atom.h"
#include "obdetails.h"

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/obiter.h>


namespace OpenBabel
{

bool AtomRoles::HasRole(const std::string &test_role) {
	return (_roles.find(test_role) != _roles.end());
}

void AtomRoles::ClearRoles() {
	_roles.clear();
}

void AtomRoles::SetRole(const std::string &possibly_new_role) {
	_roles.insert(possibly_new_role);
}

bool AtomRoles::DeleteRole(const std::string &role) {
	bool deleted = false;
	if (HasRole(role)) {
		_roles.erase(role);
	}
	return deleted;
}


Topology::Topology(OBMol *parent_mol) {
	orig_molp = parent_mol;
	simplified_net = initMOFwithUC(parent_mol);

	VirtualMol connections(&simplified_net);
	PseudoAtomMap atom_translator(&simplified_net, orig_molp);
	std::map<OBAtom*, AtomRoles> role_translator();    // initialize to empty.  Automatically will add elements

	// Initialize simplified_net via copying orig_mol and creating the 1:1 mapping
	std::map<OBAtom*, PseudoAtom> atom_to_pa;  // temporary var for bond accounting
	FOR_ATOMS_OF_MOL(orig_atom, *orig_molp) {
		OBAtom* new_atom;
		// TODO: set the simplified net to fake pseudoatoms instead of real elements
		new_atom = formAtom(&simplified_net, orig_atom->GetVector(), orig_atom->GetAtomicNum());
		atom_translator[new_atom].AddAtom(&*orig_atom);
		atom_to_pa[&*orig_atom] = new_atom;
	}
	FOR_BONDS_OF_MOL(orig_bond, *orig_molp) {
		OBAtom* orig_a1 = orig_bond->GetBeginAtom();
		OBAtom* orig_a2 = orig_bond->GetEndAtom();
		int new_order = 1;  // don't need to assign bond orders in the simplified topolgy
		formBond(&simplified_net, atom_to_pa[orig_a1], atom_to_pa[orig_a2], new_order);
	}
}

bool Topology::IsConnection(OBAtom* a) {
	return connections.HasAtom(a);
}

VirtualMol Topology::GetOrigAtomsOfRole(const std::string &role) {
	VirtualMol match(&simplified_net);
	for (std::map<PseudoAtom, AtomRoles>::iterator it=role_translator.begin(); it!=role_translator.end(); ++it) {
		if (it->second.HasRole(role)) {
			match.AddVirtualMol(atom_translator[it->first]);
		}
	}
	return match;
}

OBMol Topology::ToOBMol() {
	// TODO: eventually this code will assign pseudo-atom types based on SMILES (like ElementGen),
	// but for now, just spit out the OBMol of interest
	return simplified_net;
}

} // end namespace OpenBabel
