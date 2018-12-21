#include "topology.h"
#include "framework.h"
#include "virtual_mol.h"
#include "pseudo_atom.h"
#include "periodic.h"
#include "obdetails.h"

#include <vector>
#include <string>
#include <set>
#include <tuple>
#include <queue>

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/obiter.h>
#include <openbabel/generic.h>


namespace OpenBabel
{

bool AtomRoles::HasRole(const std::string &test_role) {
	return (_roles.find(test_role) != _roles.end());
}

void AtomRoles::ClearRoles() {
	_roles.clear();
}

void AtomRoles::AddRole(const std::string &possibly_new_role) {
	_roles.insert(possibly_new_role);
}

bool AtomRoles::RemoveRole(const std::string &role) {
	bool deleted = false;
	if (HasRole(role)) {
		_roles.erase(role);
	}
	return deleted;
}


Connections::Connections(OBMol* parent) {
	parent_net = parent;
}

void Connections::AddConn(PseudoAtom conn, PseudoAtom begin, PseudoAtom end) {
	if (begin == end) { return; }  // trivial loop: TODO raise ERROR
	std::pair<PseudoAtom, PseudoAtom> endpoints(begin, end);
	conn2endpts[conn] = endpoints;
	endpt_nbors[begin].insert(end);
	endpt_nbors[end].insert(begin);
	endpt_conns[begin].insert(conn);
	endpt_conns[end].insert(conn);
}

void Connections::RemoveConn(PseudoAtom conn) {
	std::pair<PseudoAtom, PseudoAtom> endpoints;
	endpoints = conn2endpts[conn];
	conn2endpts.erase(conn);
	endpt_nbors.erase(endpoints.first);
	endpt_nbors.erase(endpoints.second);
	endpt_conns.erase(endpoints.first);
	endpt_conns.erase(endpoints.second);
}

bool Connections::IsConn(PseudoAtom atom) {
	return (conn2endpts.find(atom) == conn2endpts.end());
}

AtomSet Connections::GetAtomConns(PseudoAtom atom) {
	return endpt_conns[atom];
}

AtomSet Connections::GetAtomNeighbors(PseudoAtom atom) {
	return endpt_nbors[atom];
}

PseudoAtom Connections::GetConn(PseudoAtom begin, PseudoAtom end) {
	// Error if the second atom is not a neighbor of the first
	AtomSet begin_nbors = GetAtomNeighbors(begin);
	if (begin_nbors.find(end) == begin_nbors.end()) { return NULL; }

	PseudoAtom connection_site = NULL;
	AtomSet all_conns = GetAtomConns(begin);
	for (AtomSet::iterator it=all_conns.begin(); it!=all_conns.end(); ++it) {
		PseudoAtom test_conn = *it;
		if ((conn2endpts[test_conn].first==end) || (conn2endpts[test_conn].second==end)) {
			connection_site = test_conn;
		}
	}
	return connection_site;
}

AtomSet Connections::GetConnEndpoints(PseudoAtom conn) {
	std::set<PseudoAtom> endpoints;
	endpoints.insert(conn2endpts[conn].first);
	endpoints.insert(conn2endpts[conn].second);
	return endpoints;
}

VirtualMol Connections::GetInternalConns(VirtualMol atoms) {
	if (atoms.GetParent() != parent_net) { return VirtualMol(NULL); };  // error
	VirtualMol int_conn(parent_net);
	AtomSet pa = atoms.GetAtoms();
	for (AtomSet::iterator a_it=pa.begin(); a_it!=pa.end(); ++a_it) {
		AtomSet pa_conns = endpt_conns[*a_it];
		for (AtomSet::iterator conn_it=pa_conns.begin(); conn_it!=pa_conns.end(); ++conn_it) {
			std::pair<PseudoAtom, PseudoAtom> endpoints = conn2endpts[*conn_it];
			if (atoms.HasAtom(endpoints.first) && atoms.HasAtom(endpoints.second)) {
				// Only add connections if both endpoints are internal
				int_conn.AddAtom(*conn_it);
			}
		}
	}
	return int_conn;
}

/*
bool Connections::CheckConsistency() {
	// This would be a nice feature eventually
	// number of bonds upstream must equal number of connections.
	// also the two internal databases should match
	return true;
}
*/


Topology::Topology(OBMol *parent_mol) {
	orig_molp = parent_mol;
	simplified_net = initMOFwithUC(parent_mol);

	Connections conns(&simplified_net);
	VirtualMol deleted_atoms(orig_molp);
	PseudoAtomMap pa_to_act(&simplified_net, orig_molp);
	std::map<OBAtom*, AtomRoles> act_roles();    // initialize to empty.  Automatically will add elements

	// Initialize simplified_net via copying orig_mol and creating the 1:1 mapping
	FOR_ATOMS_OF_MOL(orig_atom, *orig_molp) {
		OBAtom* new_atom;
		new_atom = formAtom(&simplified_net, orig_atom->GetVector(), DEFAULT_ELEMENT);
		pa_to_act[new_atom].AddAtom(&*orig_atom);
		act_to_pa[&*orig_atom] = new_atom;
	}
	// Bonds in the simplified net are handled specially with a shadow Connections table
	FOR_BONDS_OF_MOL(orig_bond, *orig_molp) {
		PseudoAtom begin_pa = act_to_pa[orig_bond->GetBeginAtom()];
		PseudoAtom end_pa = act_to_pa[orig_bond->GetEndAtom()];
		ConnectAtoms(begin_pa, end_pa);
	}
}

bool Topology::IsConnection(PseudoAtom a) {
	return conns.IsConn(a);
}

VirtualMol Topology::GetOrigAtomsOfRole(const std::string &role) {
	VirtualMol match(orig_molp);
	for (std::map<OBAtom*, AtomRoles>::iterator it=act_roles.begin(); it!=act_roles.end(); ++it) {
		if (it->second.HasRole(role)) {
			match.AddAtom(it->first);
		}
	}
	return match;
}

void Topology::SetRoleToAtoms(const std::string &role, VirtualMol atoms, bool val) {
	// Adds/removes the role from a list of original atoms in a VirtualMol
	if (atoms.GetParent() != orig_molp) {
		// Error
		return;
	}
	std::set<OBAtom*> atom_list = atoms.GetAtoms();
	for (std::set<OBAtom*>::iterator it=atom_list.begin(); it!=atom_list.end(); ++it) {
		if (val) {
			act_roles[*it].AddRole(role);
		} else {
			act_roles[*it].RemoveRole(role);
		}
	}
}

int Topology::RemoveOrigAtoms(VirtualMol atoms) {
	// FIXME: come back and implement at the end
	if (atoms.GetParent() != orig_molp) {
		return 0;  // error
	}
	std::set<OBAtom*> atom_set = atoms.GetAtoms();
	for (std::set<OBAtom*>::iterator it=atom_set.begin(); it!=atom_set.end(); ++it) {
		deleted_atoms.AddAtom(*it);
		// delete from simplified net
		// other accounting to take care of?

	}
	return 1;
	// remove from the simplified net
	//deleted_atoms
}

PseudoAtom Topology::ConnectAtoms(PseudoAtom begin, PseudoAtom end, vector3 *pos) {
	// Form the pseudo atom
	// TODO: think about if they're already connected.  That's not a problem
	// (in fact it's the whole point, for structures like MOF-5), but I should consider
	// what happens with the Connections class, etc., under those circumstances.
	// Also how we consider valence of an atom, etc.
	vector3 conn_pos;
	if (pos != NULL) {
		conn_pos = *pos;
	} else {  // use the midpoint by default
		OBUnitCell* uc = getPeriodicLattice(&simplified_net);
		vector3 begin_pos = begin->GetVector();
		vector3 end_pos = uc->UnwrapCartesianNear(end->GetVector(), begin_pos);

		conn_pos = uc->WrapCartesianCoordinate((begin_pos+end_pos) / 2.0);
	}
	PseudoAtom new_conn = formAtom(&simplified_net, conn_pos, CONNECTION_ELEMENT);

	// Form bonds and update accounting for Connections object
	formBond(&simplified_net, begin, new_conn, 1);
	formBond(&simplified_net, end, new_conn, 1);
	conns.AddConn(new_conn, begin, end);
	return new_conn;
}

void Topology::DeleteConnection(PseudoAtom conn) {
	// Removes connections between two atoms (no longer directly bonded through a connection site)
	conns.RemoveConn(conn);
	simplified_net.DeleteAtom(conn);  // automatically deletes attached bonds
}

void Topology::DeleteConnection(PseudoAtom begin, PseudoAtom end) {
	DeleteConnection(conns.GetConn(begin, end));
}

PseudoAtom Topology::CollapseOrigAtoms(VirtualMol atoms) {
	// FIXME still implementing
	if (atoms.GetParent() != orig_molp) {
		return NULL;  // error
	}

	// Find the relevant pseudoatoms and their external connections
	VirtualMol pa_int(&simplified_net);
	std::set<OBAtom*> act_atoms = atoms.GetAtoms();
	for (std::set<OBAtom*>::iterator it=act_atoms.begin(); it!=act_atoms.end(); ++it) {
		pa_int.AddAtom(act_to_pa[*it]);
	}
	// Add connections inside pa_int
	// this whole sectino will be replaced by GetInternalConns:
	std::set<OBAtom*> orig_pa = pa_int.GetAtoms();
	std::set<OBAtom*> internal_connections;
	for (std::set<OBAtom*>::iterator it=orig_pa.begin(); it!=orig_pa.end(); ++it) {
		FOR_NBORS_OF_ATOM(nbor, *it) {
			if (IsConnection(&*nbor)) {
				std::set<OBAtom*> visited;
				visited.insert(*it);
				OBAtom* next = &*nbor;
				// stuff
				// TODO: going to get connections sorted out
				//while
			}
		}
	}
	// get X_CONN - loop over conn
	// Check for consistency in the
	// loop to verify pa_int<unique atoms>.length vs. number in VirtualMol




	// TODO deleteME
	// Get external connections for the corresponding pseudo atoms
	ConnIntToExt atom_conn = atoms.GetExternalConnections();
	// warning: we can't do the external connections this way.  What if there were already simplifications to the topology?
	// However, we can probably just translate atoms -> pa, then run GetExternalConnections()
	ConnIntToExt pa_conn;
	for (ConnIntToExt::iterator it=atom_conn.begin(); it!=atom_conn.end(); ++it) {
		PseudoAtom pa_int = act_to_pa[it->first];
		PseudoAtom pa_ext = act_to_pa[it->second];
		// some sort of consistency check required
		// but what to do about connection pseudo atoms?
		std::pair<OBAtom*, OBAtom*> pa_bond(pa_int, pa_ext);
		pa_conn.insert(pa_bond);
	}
	//ConnIntToExt atom_conn = atoms.GetExternalConnections();
	// check NULL and length of pseudo atoms?
	// will also need to update the PA/act arrays

	// Put the connection point 1/3 of the way between the centroid and the connection midpoint to the exterior
	// (e.g. 1/3 of the way between a BDC centroid and the O-M bond in the -COO group).
	// In a simplified M-X-X-M' system, this will have the convenient property of being mostly equidistant.
	// Note: this follows the convention of many top-down MOF generators placing the connection point halfway on the node-linker bond.
	// In this circumstance, the convention also has the benefit that a linker with many connections to the same metal (-COO)
	// or connections to multiple metals (MOF-74 series) have unique positions for the X_CONN pseudo atoms.
	OBMol center_of_atoms = atoms.ToOBMol(false);  // TODO RENAME
	vector3 centroid = getCentroid(&center_of_atoms, false);  // without bonds or weights
	for (ConnIntToExt::iterator it=pa_conn.begin(); it!=pa_conn.end(); ++it) {
		OBUnitCell* lattice = getPeriodicLattice(&simplified_net);
		OBAtom* pa_int = it->first;
		OBAtom* pa_ext = it->second;

		vector3 int_loc = lattice->UnwrapCartesianNear(pa_int->GetVector(), centroid);
		vector3 ext_loc = lattice->UnwrapCartesianNear(pa_ext->GetVector(), int_loc);
		vector3 conn_loc = lattice->WrapCartesianCoordinate((4.0*centroid + int_loc + ext_loc) / 6.0);

		OBAtom* conn_atom = formAtom(&simplified_net, conn_loc, CONNECTION_ELEMENT);
		//connections.AddAtom(conn_atom);
		formBond(&simplified_net, conn_atom, pa_int, 1);  // Connect to internal
		formBond(&simplified_net, conn_atom, pa_ext, 1);
	}
	// update arrays, etc.


	return NULL;
}

OBMol Topology::ToOBMol() {
	// TODO: eventually this code will assign pseudo-atom types based on SMILES (like ElementGen),
	// but for now, just spit out the OBMol of interest

	// For the interim, let's try coloring the atoms as a test.
	// This will not likely be the implementation for the final version of the code, but it's worth trying now
	for (std::map<OBAtom*, AtomRoles>::iterator it=act_roles.begin(); it!=act_roles.end(); ++it) {
		PseudoAtom a = act_to_pa[it->first];
		if (it->second.HasRole("node")) {
			changeAtomElement(a, 40);  // Zr (teal)
		} else if (it->second.HasRole("linker")) {
			changeAtomElement(a, 7);  // N (blue)
		}
	}
	// I did the coloring this way out of convenience, but honestly it's actually
	// a really good way to visualize how the net turned out.
	// A web app to combine/hide the different layers could work too.
	return simplified_net;
}

} // end namespace OpenBabel
