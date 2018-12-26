#include "topology.h"
#include "framework.h"
#include "virtual_mol.h"
#include "pseudo_atom.h"
#include "periodic.h"
#include "obdetails.h"

#include <vector>
#include <string>
#include <set>
#include <utility>  // std::pair
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


ConnectionTable::ConnectionTable(OBMol* parent) {
	parent_net = parent;
}

void ConnectionTable::AddConn(PseudoAtom conn, PseudoAtom begin, PseudoAtom end) {
	if (begin == end) {
		obErrorLog.ThrowError(__FUNCTION__, "Adding a connection for a trivial loop.", obWarning);
		// TODO: degrade to an Info level notification
		// FIXME: before that, make sure the rest of the code is sufficiently robust to begin==end

		// Do not return or otherwise cause an error.  This condition can naturally occur with
		// a ConnectionTable (e.g. MOF-5) even though it can't with native OBBonds.
	}
	std::pair<PseudoAtom, PseudoAtom> endpoints(begin, end);
	conn2endpts[conn] = endpoints;
	endpt_conns[begin].insert(conn);
	endpt_conns[end].insert(conn);
}

void ConnectionTable::RemoveConn(PseudoAtom conn) {
	std::pair<PseudoAtom, PseudoAtom> endpoints;
	endpoints = conn2endpts[conn];
	conn2endpts.erase(conn);
	endpt_conns[endpoints.first].erase(conn);
	endpt_conns[endpoints.second].erase(conn);
}

bool ConnectionTable::IsConn(PseudoAtom atom) {
	return (conn2endpts.find(atom) != conn2endpts.end());
}

AtomSet ConnectionTable::GetAtomConns(PseudoAtom endpt) {
	return endpt_conns[endpt];
}

bool ConnectionTable::HasNeighbor(PseudoAtom begin, PseudoAtom end) {
	AtomSet begin_conn = GetAtomConns(begin);
	for (AtomSet::iterator it=begin_conn.begin(); it!=begin_conn.end(); ++it) {
		if ((conn2endpts[*it].first == end) || (conn2endpts[*it].second == end)) {
			return true;
		}
	}
	return false;
}

AtomSet ConnectionTable::GetConnEndpointSet(PseudoAtom conn) {
	std::set<PseudoAtom> endpoints;
	endpoints.insert(conn2endpts[conn].first);
	endpoints.insert(conn2endpts[conn].second);
	return endpoints;
}

std::pair<PseudoAtom, PseudoAtom> ConnectionTable::GetConnEndpoints(PseudoAtom conn) {
	// unordered pair of <begin, end> endpoints
	return conn2endpts[conn];
}

VirtualMol ConnectionTable::GetInternalConns(VirtualMol atoms) {
	if (atoms.GetParent() != parent_net) {
		obErrorLog.ThrowError(__FUNCTION__, "VirtualMol parent mismatch", obWarning);
		return VirtualMol();  // null parent
	};
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
bool ConnectionTable::CheckConsistency() {
	// This would be a nice feature eventually
	// number of bonds upstream must equal number of connections.
	// also the two internal databases should match
	return true;
}
*/


Topology::Topology(OBMol *parent_mol) {
	orig_molp = parent_mol;
	simplified_net = initMOFwithUC(parent_mol);

	// Remember not to declare the object types in the constructor.
	// We're trying to initialize class members, not declare local variables with the same name.
	conns = ConnectionTable(&simplified_net);
	deleted_atoms = VirtualMol(orig_molp);
	pa_to_act = PseudoAtomMap(&simplified_net, orig_molp);
	pa_roles = std::map<OBAtom*, AtomRoles>();    // initialize to empty.  Automatically will add elements

	// Initialize simplified_net via copying orig_mol and creating the 1:1 mapping
	FOR_ATOMS_OF_MOL(orig_atom, *orig_molp) {
		OBAtom* new_atom;
		new_atom = formAtom(&simplified_net, orig_atom->GetVector(), DEFAULT_ELEMENT);
		pa_to_act[new_atom] = VirtualMol(&*orig_atom);
		act_to_pa[&*orig_atom] = new_atom;
	}
	// Bonds in the simplified net are handled specially with a shadow ConnectionTable object
	FOR_BONDS_OF_MOL(orig_bond, *orig_molp) {
		PseudoAtom begin_pa = act_to_pa[orig_bond->GetBeginAtom()];
		PseudoAtom end_pa = act_to_pa[orig_bond->GetEndAtom()];
		ConnectAtoms(begin_pa, end_pa);
	}
}

bool Topology::IsConnection(PseudoAtom a) {
	return conns.IsConn(a);
}

VirtualMol Topology::GetAtomsOfRole(const std::string &role) {
	VirtualMol match(&simplified_net);
	for (std::map<OBAtom*, AtomRoles>::iterator it=pa_roles.begin(); it!=pa_roles.end(); ++it) {
		if (it->second.HasRole(role)) {
			match.AddAtom(it->first);
		}
	}
	return match;
}

VirtualMol Topology::GetAtoms(bool include_conn) {
	VirtualMol atoms(&simplified_net);
	FOR_ATOMS_OF_MOL(a, simplified_net) {
		if (include_conn || !IsConnection(&*a)) {
			atoms.AddAtom(&*a);
		}
	}
	return atoms;
}

bool Topology::AtomHasRole(PseudoAtom atom, const std::string &role) {
	if (pa_roles.find(atom) == pa_roles.end()) {
		obErrorLog.ThrowError(__FUNCTION__, "Could not find pseudoatom in pa_roles mapping", obWarning);
		return false;
	}
	return pa_roles[atom].HasRole(role);
}

void Topology::SetRoleToAtom(const std::string &role, PseudoAtom atom, bool val) {
	if (val) {
		pa_roles[atom].AddRole(role);
	} else {
		pa_roles[atom].RemoveRole(role);
	}
}

void Topology::SetRoleToAtoms(const std::string &role, VirtualMol atoms, bool val) {
	// Adds/removes the role from a list of PseudoAtoms in the simplified net
	if (atoms.GetParent() != &simplified_net) {
		obErrorLog.ThrowError(__FUNCTION__, "VirtualMol needs to contain PseudoAtoms of the simplified net.", obError);
		return;
	}
	std::set<OBAtom*> atom_list = atoms.GetAtoms();
	for (std::set<OBAtom*>::iterator it=atom_list.begin(); it!=atom_list.end(); ++it) {
		SetRoleToAtom(role, *it, val);
	}
}

int Topology::RemoveOrigAtoms(VirtualMol atoms) {
	// FIXME: come back and implement at the end
	// Used for solvent removal, etc.
	if (atoms.GetParent() != orig_molp) {
		//obErrorLog.ThrowError(__FUNCTION__, "VirtualMol needs to contain child atoms of the original, unsimplified MOF", obError);
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

VirtualMol Topology::OrigToPseudo(VirtualMol orig_atoms) {
	if (orig_atoms.GetParent() != orig_molp) {
		obErrorLog.ThrowError(__FUNCTION__, "VirtualMol needs to contain child atoms of the original, unsimplified MOF", obError);
		return VirtualMol();
	}
	std::set<OBAtom*> act_atoms = orig_atoms.GetAtoms();
	VirtualMol pa(&simplified_net);

	// Find the relevant set of pseudoatoms
	for (std::set<OBAtom*>::iterator it=act_atoms.begin(); it!=act_atoms.end(); ++it) {
		pa.AddAtom(act_to_pa[*it]);
	}
	// TODO consider a consistency check that the PA's don't include any other atoms (a length check for fragment vs. sum of PA AtomSets)

	return pa;
}

VirtualMol Topology::PseudoToOrig(VirtualMol pa_atoms) {
	if (pa_atoms.GetParent() != &simplified_net) {
		obErrorLog.ThrowError(__FUNCTION__, "VirtualMol needs to contain child atoms of the simplified net.", obError);
		return VirtualMol();
	}

	VirtualMol orig_atoms(orig_molp);
	AtomSet pa_set = pa_atoms.GetAtoms();
	for (std::set<OBAtom*>::iterator it=pa_set.begin(); it!=pa_set.end(); ++it) {
		orig_atoms.AddVirtualMol(pa_to_act[*it]);
	}
	return orig_atoms;
}

PseudoAtom Topology::ConnectAtoms(PseudoAtom begin, PseudoAtom end, vector3 *pos) {
	// Form the pseudo atom
	// It's not a problem if they're already connected.
	// In fact it's the whole point for structures like MOF-5, and the ConnectionTable is built to handle it.
	// As such, the valence of a pseudoatom will be it's OBAtom valence or size of ConnectionTable::GetAtomConns
	vector3 conn_pos;
	if (pos != NULL) {
		conn_pos = *pos;
	} else {  // use the midpoint by default
		if (begin == end) {
			obErrorLog.ThrowError(__FUNCTION__, "Underdefined: need an explicit connection location to connect an atom to itself.", obError);
			return NULL;
		}

		OBUnitCell* uc = getPeriodicLattice(&simplified_net);
		vector3 begin_pos = begin->GetVector();
		vector3 end_pos = uc->UnwrapCartesianNear(end->GetVector(), begin_pos);

		conn_pos = uc->WrapCartesianCoordinate((begin_pos+end_pos) / 2.0);
	}
	PseudoAtom new_conn = formAtom(&simplified_net, conn_pos, CONNECTION_ELEMENT);

	// Form bonds and update accounting
	formBond(&simplified_net, begin, new_conn, 1);
	formBond(&simplified_net, end, new_conn, 1);
	conns.AddConn(new_conn, begin, end);
	pa_roles[new_conn].AddRole("connection");
	pa_to_act[new_conn] = VirtualMol(orig_molp);

	return new_conn;
}

void Topology::DeleteConnection(PseudoAtom conn) {
	// Removes connections between two atoms (no longer directly bonded through a connection site)
	conns.RemoveConn(conn);
	simplified_net.DeleteAtom(conn);  // automatically deletes attached bonds
	pa_roles.erase(conn);
	pa_to_act.RemoveAtom(conn);
}

void Topology::DeleteAtomAndConns(PseudoAtom atom) {
	// Removes an atom and relevant bonds/connections
	AtomSet nbors;
	FOR_NBORS_OF_ATOM(nbor, *atom) {
		if (!IsConnection(&*nbor)) {
			obErrorLog.ThrowError(__FUNCTION__, "Undefined behavior: Found atoms directly bonded to each other (instead of an intermediate connection site).", obError);
			return;
		}
		nbors.insert(&*nbor);
	}
	for (AtomSet::iterator it=nbors.begin(); it!=nbors.end(); ++it) {
		DeleteConnection(*it);
	}
	simplified_net.DeleteAtom(atom);  // automatically deletes bonds

	// Remove original atoms if present
	AtomSet act_atoms = pa_to_act[atom].GetAtoms();
	for (AtomSet::iterator it=act_atoms.begin(); it!=act_atoms.end(); ++it) {
		act_to_pa[*it] = NULL;
		deleted_atoms.AddAtom(*it);
	}
	pa_to_act.RemoveAtom(atom);  // and remove it from the key of PA's
	pa_roles.erase(atom);
}

ConnIntToExt Topology::GetConnectedAtoms(VirtualMol internal_pa) {
	// Gets the next shell of PseudoAtom neighbors external to the internal VirtualMol.
	// Automatically passes over connection pseudoatoms.

	// VirtualMol::GetExternalBonds() will also return internal connections, so let's specify those ahead of time
	VirtualMol pa_with_conns = FragmentWithIntConns(internal_pa);
	// Then GetExternalBonds() can only return extra bonds
	ConnIntToExt bonds = pa_with_conns.GetExternalBonds();

	// Translate external connections to external pseudoatoms
	ConnIntToExt external_nbors;
	for (ConnIntToExt::iterator e_it=bonds.begin(); e_it!=bonds.end(); ++e_it) {
		PseudoAtom int_to_conn = e_it->first;
		PseudoAtom ext_conn = e_it->second;
		PseudoAtom ext_from_conn = NULL;

		// Figure out which of the conn_endpts is the new external endpoint
		std::pair<PseudoAtom, PseudoAtom> conn_endpts = conns.GetConnEndpoints(ext_conn);
		if (int_to_conn == conn_endpts.first) {
			ext_from_conn = conn_endpts.second;
		} else if (int_to_conn == conn_endpts.second) {
			ext_from_conn = conn_endpts.first;
		} else {
			obErrorLog.ThrowError(__FUNCTION__, "AssertionError: ended at an unexpected endpoint", obError);
			return ConnIntToExt();
		}
		std::pair<OBAtom*, OBAtom*> ConnExt(int_to_conn, ext_from_conn);
		external_nbors.insert(ConnExt);
	}

	return external_nbors;
}

PseudoAtom Topology::CollapseFragment(VirtualMol pa_fragment) {
	// Simplifies the net by combining pa_fragment PseudoAtoms into a single point,
	// maintaining existing connections.  Returns a pointer to the generated PseudoAtom.
	if (pa_fragment.GetParent() != &simplified_net) {
		obErrorLog.ThrowError(__FUNCTION__, "VirtualMol needs to contain child atoms of the simplified net.", obError);
		return NULL;
	}

	// Get external neighbors on the other side of the connections.
	ConnIntToExt external_nbors = GetConnectedAtoms(pa_fragment);

	// Make the new pseudoatom at the centroid
	OBMol mol_orig_pa = FragmentToOBMolNoConn(pa_fragment);  // atoms and bonds
	vector3 centroid = getCentroid(&mol_orig_pa, false);  // without weights
	PseudoAtom new_atom = formAtom(&simplified_net, centroid, DEFAULT_ELEMENT);

	// Put the connection point 1/3 of the way between the centroid and the connection midpoint to the exterior
	// (e.g. 1/3 of the way between a BDC centroid and the O-M bond in the -COO group).
	// In a simplified M-X-X-M' system, this will have the convenient property of being mostly equidistant.
	// Note: this follows the convention of many top-down MOF generators placing the connection point halfway on the node-linker bond.
	// In this circumstance, the convention also has the benefit that a linker with many connections to the same metal (-COO)
	// or connections to multiple metals (MOF-74 series) have unique positions for the X_CONN pseudo atoms.
	for (ConnIntToExt::iterator it=external_nbors.begin(); it!=external_nbors.end(); ++it) {
		OBAtom* pa_int = it->first;
		OBAtom* pa_ext = it->second;

		// Positions of the internal/external bonding atoms, to get the bond location
		OBUnitCell* lattice = getPeriodicLattice(&simplified_net);
		vector3 int_loc = lattice->UnwrapCartesianNear(pa_int->GetVector(), centroid);
		vector3 ext_loc = lattice->UnwrapCartesianNear(pa_ext->GetVector(), int_loc);

		vector3 conn_loc = lattice->WrapCartesianCoordinate((4.0*centroid + int_loc + ext_loc) / 6.0);
		ConnectAtoms(new_atom, pa_ext, &conn_loc);
	}

	// Update the mapping between PA's and original atoms, then delete the
	// original pseudoatoms corresponding with the fragment
	VirtualMol pa_without_conn = FragmentWithoutConns(pa_fragment);
	AtomSet act_atoms = PseudoToOrig(pa_without_conn).GetAtoms();
	for (AtomSet::iterator it=act_atoms.begin(); it!=act_atoms.end(); ++it) {
		act_to_pa[*it] = new_atom;
	}
	AtomSet orig_pa_set = pa_without_conn.GetAtoms();
	for (AtomSet::iterator it=orig_pa_set.begin(); it!=orig_pa_set.end(); ++it) {
		pa_to_act.RemoveAtom(*it);
		DeleteAtomAndConns(*it);
	}

	return new_atom;
}

void Topology::MergeAtomToAnother(PseudoAtom from, PseudoAtom to) {
	// like CollapseFragment, but moves contents from a "from"
	// PseudoAtom to a compatible "to" atom without moving the second atom.

	// Check compatibility of the two atoms:
	// all of the first PA's neighbors must be either the to/from atom or a mutual neighbor.
	FOR_NBORS_OF_ATOM(from_conn, *from) {
		FOR_NBORS_OF_ATOM(from_nbor, *from_conn) {
			if ((&*from_nbor!=to) && (&*from_nbor!=from) && !conns.HasNeighbor(to, &*from_nbor)) {
				obErrorLog.ThrowError(__FUNCTION__, "Cannot merge incompatible atoms with mismatched neighbors", obError);
				return;
			}
		}
	}

	// Move content from origin to destination, then delete the origin
	AtomSet from_mapping = pa_to_act[from].GetAtoms();
	for (AtomSet::iterator it=from_mapping.begin(); it!=from_mapping.end(); ++it) {
		if (pa_to_act[to].HasAtom(*it)) {
			obErrorLog.ThrowError(__FUNCTION__, "AssertionError: orig_atom should not be a child of both origin and destination PseudoAtoms", obError);
		}
		pa_to_act[to].AddAtom(*it);

		if (act_to_pa[*it] != from) {
			obErrorLog.ThrowError(__FUNCTION__, "AssertionError: inconsistency in original PA ownership of orig_atom", obError);
		}
		act_to_pa[*it] = to;
	}
	// skip updating the atom roles
	DeleteAtomAndConns(from);
}

OBMol Topology::FragmentToOBMolNoConn(VirtualMol pa_fragment) {
	// Based on VirtualMol::ToOBMol
	OBMol mol = initMOFwithUC(&simplified_net);
	std::map<OBAtom*, OBAtom*> virtual_to_mol;
	AtomSet fragment_atoms = pa_fragment.GetAtoms();
	// Copy atoms
	for (AtomSet::iterator it=fragment_atoms.begin(); it!=fragment_atoms.end(); ++it) {
		OBAtom* fragment_atom = (*it);
		if (IsConnection(fragment_atom)) { continue; }  // skip over connections
		OBAtom* mol_atom = formAtom(&mol, fragment_atom->GetVector(), fragment_atom->GetAtomicNum());
		virtual_to_mol[fragment_atom] = mol_atom;
	}

	// Convert connections into bonds
	AtomSet internal_conns = conns.GetInternalConns(pa_fragment).GetAtoms();
	for (AtomSet::iterator it=internal_conns.begin(); it!=internal_conns.end(); ++it) {
		std::pair<PseudoAtom, PseudoAtom> begin_end = conns.GetConnEndpoints(*it);
		PseudoAtom begin = virtual_to_mol[begin_end.first];
		PseudoAtom end = virtual_to_mol[begin_end.second];
		if (mol.GetBond(begin, end)) {
			obErrorLog.ThrowError(__FUNCTION__, "OBMol will have undefined connections due to multiply-bonded begin/end atoms.", obWarning);
		} else {
			formBond(&mol, begin, end, 1);
		}
	}

	return mol;
}

OBMol Topology::ToOBMol() {
	// TODO: eventually this code will assign pseudo-atom types based on SMILES (like ElementGen),
	// but for now, just spit out the OBMol of interest.
	// In the final version, be sure to initMOFwithUC to get the lattice params correctly

	// For the interim, let's try coloring the atoms as a test.
	// This will not likely be the implementation for the final version of the code, but it's worth trying now
	for (std::map<OBAtom*, AtomRoles>::iterator it=pa_roles.begin(); it!=pa_roles.end(); ++it) {
		PseudoAtom a = act_to_pa[it->first];
		if (it->second.HasRole("node")) {
			changeAtomElement(it->first, 40);  // Zr (teal)
		} else if (it->second.HasRole("linker")) {
			changeAtomElement(it->first, 7);  // N (blue)
		} else if (it->second.HasRole("connection")) {
			changeAtomElement(it->first, 8);  // O (red)
		}
	}
	// I did the coloring this way out of convenience, but honestly it's actually
	// a really good way to visualize how the net turned out.
	// A web app to combine/hide the different layers could work too.
	return simplified_net;
}

VirtualMol Topology::FragmentWithoutConns(VirtualMol fragment) {
	// remove connection pseudoatoms from a VirtualMol
	VirtualMol cleaned(fragment.GetParent());
	AtomSet all_atoms = fragment.GetAtoms();
	for (AtomSet::iterator it=all_atoms.begin(); it!=all_atoms.end(); ++it) {
		if (!IsConnection(*it)) {
			cleaned.AddAtom(*it);
		}
	}
	return cleaned;
}

VirtualMol Topology::FragmentWithIntConns(VirtualMol fragment) {
	VirtualMol pa_and_conns = fragment;
	pa_and_conns.AddVirtualMol(conns.GetInternalConns(fragment));
	return pa_and_conns;
}

} // end namespace OpenBabel
