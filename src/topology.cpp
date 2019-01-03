#include "topology.h"
#include "framework.h"
#include "virtual_mol.h"
#include "pseudo_atom.h"
#include "periodic.h"
#include "obdetails.h"
#include "invector.h"

#include <iostream>
#include <fstream>
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

PseudoAtom ConnectionTable::GetOtherEndpoint(PseudoAtom conn, PseudoAtom begin) {
	// Given a connection and one endpoint, what is the other endpoint?
	std::pair<PseudoAtom, PseudoAtom> endpoints = GetConnEndpoints(conn);
	if (begin == endpoints.first) {
		return endpoints.second;
	} else if (begin == endpoints.second) {
		return endpoints.first;
	} else {
		obErrorLog.ThrowError(__FUNCTION__, "PseudoAtom not part of the connection", obError);
		return NULL;
	}
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
	deleted_atoms = std::map<std::string, VirtualMol>();  // empty: initially all atoms from orig_mol exist
	pa_to_act = PseudoAtomMap(&simplified_net, orig_molp);
	pa_roles = std::map<OBAtom*, std::string>();    // initialize to empty.  Automatically will add elements

	// Initialize simplified_net via copying orig_mol and creating the 1:1 mapping
	FOR_ATOMS_OF_MOL(orig_atom, *orig_molp) {
		OBAtom* new_atom;
		new_atom = formAtom(&simplified_net, orig_atom->GetVector(), DEFAULT_ELEMENT);
		pa_to_act[new_atom] = VirtualMol(&*orig_atom);
		act_to_pa[&*orig_atom] = new_atom;
		pa_roles[new_atom] = "Original copy";
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

VirtualMol Topology::GetDeletedOrigAtoms(const std::string &deletion_reason) {
	if (deletion_reason != ALL_DELETED_ORIG_ATOMS) {
		if (deleted_atoms.find(deletion_reason) == deleted_atoms.end()) {
			obErrorLog.ThrowError(__FUNCTION__, "No deleted atoms of type: " + deletion_reason, obInfo);
			return VirtualMol(orig_molp);  // avoid [], which would make a new molecule
		} else {
			return deleted_atoms[deletion_reason];
		}
	} else {
		VirtualMol combined_deleted_atoms(orig_molp);
		for (std::map<std::string,VirtualMol>::iterator it=deleted_atoms.begin(); it!=deleted_atoms.end(); ++it) {
			combined_deleted_atoms.AddVirtualMol(it->second);
		}
		return combined_deleted_atoms;
	}
}

VirtualMol Topology::GetAtomsOfRole(const std::string &role) {
	VirtualMol match(&simplified_net);
	for (std::map<OBAtom*, std::string>::iterator it=pa_roles.begin(); it!=pa_roles.end(); ++it) {
		if (it->second == role) {
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

VirtualMol Topology::GetConnectors() {
	VirtualMol atoms(&simplified_net);
	FOR_ATOMS_OF_MOL(a, simplified_net) {
		if (IsConnection(&*a)) {
			atoms.AddAtom(&*a);
		}
	}
	return atoms;
}

bool Topology::AtomHasRole(PseudoAtom atom, const std::string &role) {
	return (GetRoleFromAtom(atom) == role);
}

void Topology::SetRoleToAtom(const std::string &role, PseudoAtom atom) {
	pa_roles[atom] = role;
}

void Topology::SetRoleToAtoms(const std::string &role, VirtualMol atoms) {
	// Adds/removes the role from a list of PseudoAtoms in the simplified net
	if (atoms.GetParent() != &simplified_net) {
		obErrorLog.ThrowError(__FUNCTION__, "VirtualMol needs to contain PseudoAtoms of the simplified net.", obError);
		return;
	}
	std::set<OBAtom*> atom_list = atoms.GetAtoms();
	for (std::set<OBAtom*>::iterator it=atom_list.begin(); it!=atom_list.end(); ++it) {
		SetRoleToAtom(role, *it);
	}
}

std::string Topology::GetRoleFromAtom(PseudoAtom atom) {
	if (pa_roles.find(atom) == pa_roles.end()) {
		obErrorLog.ThrowError(__FUNCTION__, "Unknown atom identity", obWarning);
		return "";
	}
	return pa_roles[atom];
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
	pa_roles[new_conn] = "connection";
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

void Topology::DeleteAtomAndConns(PseudoAtom atom, const std::string &role_for_orig_atoms) {
	// Removes an atom and relevant bonds/connections.
	// Also handles accounting for the orig_mol.  By default (role_for_orig_atoms="obError"), raise
	// an error if the PA still has any original atoms assigned to it.  Otherwise, assign the orig
	// atoms to a deleted_atom having the type "role_for_orig_atoms"
	if (IsConnection(atom)) {
		obErrorLog.ThrowError(__FUNCTION__, "Unexpectedly trying to delete a connection site.  Skipping deletion.", obWarning);
		return;
	}

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
	if (act_atoms.size()) {
		if (role_for_orig_atoms == DELETE_ORIG_ATOM_ERROR) {
			obErrorLog.ThrowError(__FUNCTION__, "Unexpectedly deleting a PA containing original atoms.  Assigning them to deleted_atoms[\"" + DELETE_ORIG_ATOM_ERROR + "\"]", obWarning);
		}
		if (deleted_atoms.find(role_for_orig_atoms) == deleted_atoms.end()) {
			deleted_atoms[role_for_orig_atoms] = VirtualMol(orig_molp);  // initialize if new
		}
		for (AtomSet::iterator it=act_atoms.begin(); it!=act_atoms.end(); ++it) {
			act_to_pa[*it] = NULL;
			deleted_atoms[role_for_orig_atoms].AddAtom(*it);
		}
	}
	pa_to_act.RemoveAtom(atom);  // and remove it from the key of PA's
	pa_roles.erase(atom);
}

PseudoAtom Topology::CollapseFragment(VirtualMol pa_fragment) {
	// Simplifies the net by combining pa_fragment PseudoAtoms into a single point,
	// maintaining existing connections.  Returns a pointer to the generated PseudoAtom.
	if (pa_fragment.GetParent() != &simplified_net) {
		obErrorLog.ThrowError(__FUNCTION__, "VirtualMol needs to contain child atoms of the simplified net.", obError);
		return NULL;
	}

	// Get external connection sites to the fragment
	// VirtualMol::GetExternalBondsOrConns() will also return internal connections,
	// so let's specify those ahead of time.
	VirtualMol pa_with_conns = FragmentWithIntConns(pa_fragment);
	// Then GetExternalBondsOrConns() can only return extra bonds
	ConnIntToExt external_conns = pa_with_conns.GetExternalBondsOrConns();

	// Make the new pseudoatom at the centroid
	OBMol mol_orig_pa = FragmentToOBMolNoConn(pa_fragment);  // atoms and bonds
	vector3 centroid = getCentroid(&mol_orig_pa, false);  // without weights
	PseudoAtom new_atom = formAtom(&simplified_net, centroid, DEFAULT_ELEMENT);

	// Put the connection point 1/2 way between the centroid and external connection PA
	for (ConnIntToExt::iterator it=external_conns.begin(); it!=external_conns.end(); ++it) {
		OBAtom* pa_int = it->first;
		OBAtom* pa_conn = it->second;
		OBAtom* pa_ext = conns.GetOtherEndpoint(pa_conn, pa_int);

		// Position of the new bond relative to the centroid and older, original location of the connection
		OBUnitCell* lattice = getPeriodicLattice(&simplified_net);
		vector3 old_conn_loc = lattice->UnwrapCartesianNear(pa_conn->GetVector(), centroid);
		vector3 conn_loc = lattice->WrapCartesianCoordinate((centroid + old_conn_loc) / 2.0);

		ConnectAtoms(new_atom, pa_ext, &conn_loc);
	}

	// Update the mapping between PA's and original atoms, then delete the
	// original pseudoatoms corresponding with the fragment
	VirtualMol pa_without_conn = FragmentWithoutConns(pa_fragment);
	AtomSet act_atoms = PseudoToOrig(pa_without_conn).GetAtoms();
	pa_to_act[new_atom] = VirtualMol(orig_molp);
	for (AtomSet::iterator it=act_atoms.begin(); it!=act_atoms.end(); ++it) {
		act_to_pa[*it] = new_atom;
		pa_to_act[new_atom].AddAtom(*it);
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
		pa_to_act[from].RemoveAtom(*it);

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
	for (std::map<OBAtom*, std::string>::iterator it=pa_roles.begin(); it!=pa_roles.end(); ++it) {
		PseudoAtom a = act_to_pa[it->first];
		if (it->second == "node") {
			changeAtomElement(it->first, 40);  // Zr (teal)
		} else if (it->second == "linker") {
			changeAtomElement(it->first, 7);  // N (blue)
		} else if (it->second == "connection") {
			changeAtomElement(it->first, 8);  // O (red)
		}
	}
	// I did the coloring this way out of convenience, but honestly it's actually
	// a really good way to visualize how the net turned out.
	// A web app to combine/hide the different layers could work too.
	return simplified_net;
}

void Topology::WriteSystre(const std::string &filepath, bool write_centers, bool simplify_two_conn) {
	// Write the simplified molecule to Systre for topological determination.
	// By default, this routine will account for two-connected vertices in the graph.
	// Can also print the (optional) edge_center field.
	std::ofstream ofs;
	ofs.open(filepath.c_str());

	// Write header for the molecule
	OBUnitCell* uc = getPeriodicLattice(&simplified_net);
	std::string indent = "  ";
	ofs << "# CGD file generated by Open Babel " << BABEL_VERSION << ", see http://openbabel.sf.net" << std::endl;
	ofs << "CRYSTAL" << std::endl;
	ofs << indent << "NAME " << simplified_net.GetTitle() << std::endl;
	ofs << indent << "GROUP P1" << std::endl;  // TODO: verify that Open Babel converts everything to P1
	ofs << indent << "CELL "
		<< uc->GetA() << " " << uc->GetB() << " " << uc->GetC() << " "
		<< uc->GetAlpha() << " " << uc->GetBeta() << " " << uc->GetGamma()
		<< std::endl;

	// Handle two-coordinated nodes
	VirtualMol two_coordinated(&simplified_net);
	VirtualMol multi_coordinated = GetAtoms(false);  // start with all the PA's and remove them
	VirtualMol two_xs(&simplified_net);  // connections to 2-coordinated PA's
	VirtualMol multi_xs = GetConnectors();

	if (simplify_two_conn) {
		FOR_ATOMS_OF_MOL(a, simplified_net) {
			if (IsConnection(&*a)) { continue; }
			if (a->GetValence() == 2) {
				two_coordinated.AddAtom(&*a);
				multi_coordinated.RemoveAtom(&*a);
				AtomSet connectors = conns.GetAtomConns(&*a);
				for (AtomSet::iterator it=connectors.begin(); it!=connectors.end(); ++it) {
					two_xs.AddAtom(*it);
					multi_xs.RemoveAtom(*it);
				}
			}
		}
	}  // else (if !simplify_two_conn), then two_coordinated will be empty

	int current_node = 0;
	VirtualMol visited_conns(&simplified_net);
	AtomSet multi_coordinated_set = multi_coordinated.GetAtoms();
	for (AtomSet::iterator node=multi_coordinated_set.begin(); node!=multi_coordinated_set.end(); ++node) {
		++current_node;
		vector3 frac_coords = uc->CartesianToFractional((*node)->GetVector());
		ofs << indent << "NODE " << current_node
			<< " " << (*node)->GetValence()  // coordination of the atom
			<< " " << frac_coords[0]
			<< " " << frac_coords[1]
			<< " " << frac_coords[2]
			<< std::endl;
	}

	// Export graph edges, excluding two-coordinated vertices
	std::stringstream edge_centers;
	AtomSet multi_xs_set = multi_xs.GetAtoms();
	for (AtomSet::iterator x_it=multi_xs_set.begin(); x_it!=multi_xs_set.end(); ++x_it) {
		// iterating over connectors instead of bonds
		PseudoAtom x = *x_it;
		PseudoAtom a = conns.GetConnEndpoints(x).first;
		PseudoAtom b = conns.GetConnEndpoints(x).second;

		vector3 pos_a = uc->CartesianToFractional(a->GetVector());
		vector3 pos_x = uc->UnwrapFractionalNear(uc->CartesianToFractional(x->GetVector()), pos_a);
		vector3 pos_b = uc->UnwrapFractionalNear(uc->CartesianToFractional(b->GetVector()), pos_x);
		ofs << indent << "EDGE  "
			<< pos_a[0] << " " << pos_a[1] << " " << pos_a[2] << "   "
			<< pos_b[0] << " " << pos_b[1] << " " << pos_b[2] << std::endl;
		edge_centers << indent << "# EDGE_CENTER  "
			<< pos_x[0] << " " << pos_x[1] << " " << pos_x[2] << std::endl;
	}

	// Translate two-coordinated vertices into edges
	AtomSet two_set = two_coordinated.GetAtoms();
	for (AtomSet::iterator c2_it=two_set.begin(); c2_it!=two_set.end(); ++c2_it) {
		PseudoAtom c2_linker = *c2_it;
		vector3 c2_pos = uc->CartesianToFractional(c2_linker->GetVector());
		std::vector<vector3> v2_pos;
		FOR_NBORS_OF_ATOM(c2x, *c2_linker) {
			vector3 x_pos = uc->UnwrapFractionalNear(uc->CartesianToFractional(c2x->GetVector()), c2_pos);
			vector3 vertex_pos = conns.GetOtherEndpoint(&*c2x, c2_linker)->GetVector();
			vertex_pos = uc->UnwrapFractionalNear(uc->CartesianToFractional(vertex_pos), x_pos);
			v2_pos.push_back(vertex_pos);
		}
		ofs << indent << "EDGE  "
			<< v2_pos[0][0] << " " << v2_pos[0][1] << " " << v2_pos[0][2] << "   "
			<< v2_pos[1][0] << " " << v2_pos[1][1] << " " << v2_pos[1][2] << std::endl;
		edge_centers << indent << "# EDGE_CENTER  "
			<< c2_pos[0] << " " << c2_pos[1] << " " << c2_pos[2] << std::endl;
	}

	if (write_centers) {
		ofs << edge_centers.str();
	}

	ofs << "END" << std::endl;
	ofs.close();
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

int Topology::SimplifyAxB() {
	// Remove pairs of redundant connections x1 and x2 linking together the same
	// pseudoatoms A and B (e.g. nodes and linkers) in the same direction,
	// simplifying them with a new connection site x3 at their midpoint.
	// Returns the number of modifications to connection sites.
	// Based on simplifyLX in the previous version of the code.
	// May require multiple passes to fully simplify the network (when it returns 0).

	std::vector<PseudoAtom> to_delete;  // X's to delete at the end

	AtomSet a_atoms = GetAtoms(false).GetAtoms();  // get non-connector atoms
	for (AtomSet::iterator a_it=a_atoms.begin(); a_it!=a_atoms.end(); ++a_it) {
		PseudoAtom a = *a_it;  // looping over A sites
		std::map<PseudoAtom, AtomSet> nbor_to_xs;  // all the connection X's per nbor
		AtomSet a_x_list = conns.GetAtomConns(a);
		for (AtomSet::iterator x=a_x_list.begin(); x!=a_x_list.end(); ++x) {
			PseudoAtom external_nbor = conns.GetOtherEndpoint(*x, a);
			nbor_to_xs[external_nbor].insert(*x);
		}

		for (std::map<PseudoAtom, AtomSet>::iterator b_it=nbor_to_xs.begin(); b_it!=nbor_to_xs.end(); ++b_it) {
			PseudoAtom b = b_it->first;  // looping over neighbor B sites
			AtomSet ab_xs = b_it->second;  // connections for A-x-B
			if (ab_xs.size() == 1) continue;  // no duplicate X's to check

			for (AtomSet::iterator x1=ab_xs.begin(); x1!=ab_xs.end(); ++x1) {
				for (AtomSet::iterator x2=x1; x2!=ab_xs.end(); ++x2) {
					if (
						*x1 != *x2 &&
						!inVector<PseudoAtom>(*x1, to_delete) &&
						!inVector<PseudoAtom>(*x2, to_delete)
					) {
						// Form a test molecule with A-x1-B-x2-A'
						VirtualMol test_xs(a->GetParent());
						test_xs.AddAtom(a);
						test_xs.AddAtom(b);
						test_xs.AddAtom(*x1);
						test_xs.AddAtom(*x2);
						OBMol test_xs_mol = test_xs.ToOBMol();
						if (!isPeriodicChain(&test_xs_mol)) {
							// If test_xs is periodic, then A' is in a different UC than A,
							// so X1 and X2 are a bridge.  If non-periodic (this case),
							// then X1 and X2 are redundant connections between A and B.
							to_delete.push_back(*x1);
							to_delete.push_back(*x2);
							vector3 loc = getMidpoint(*x1, *x2, false);
							ConnectAtoms(a, b, &loc);
						}
					}
				}
			}
		}
	}

	// Delete the redundant X's (which will also remove their X-L and X-M bonds)
	for (std::vector<PseudoAtom>::iterator it=to_delete.begin(); it!=to_delete.end(); ++it) {
		DeleteConnection(*it);
	}

	return to_delete.size();
}

int Topology::SplitFourVertexIntoTwoThree(PseudoAtom site) {
	// Transforms four-connected atoms to two three-connected pseudo atoms.
	// This step satisfies the convention used for MIL-47-like topologies.
	// All of the original atoms will be randomly assigned to one of the two PA's.
	// This function returns the number of simplified linkers
	if (site->GetValence() != 4) {
		obErrorLog.ThrowError(__FUNCTION__, "Cannot process site with incorrect valence.", obError);
		return 0;
	}

	// Collect the relevant pseudoatoms for the four connection sites.
	// These vectors are in sync, but in no particular order overall,
	// so vector[0] will be one of the four connections pseudorandomly.
	std::vector<PseudoAtom> site_conns;
	std::map<PseudoAtom, PseudoAtom> site_x_to_nbors;
	FOR_NBORS_OF_ATOM(nbor, *site) {
		site_conns.push_back(&*nbor);
		site_x_to_nbors[&*nbor] = conns.GetOtherEndpoint(&*nbor, site);
	}

	// The four-connected node will have two sides to break apart the one node into two.
	std::vector<PseudoAtom> side1;  // connection sites on the first site
	side1.push_back(site_conns[0]);
	side1.push_back(minAngleNbor(site, side1[0]));
	std::vector<PseudoAtom> side2;
	for (std::vector<PseudoAtom>::iterator it=site_conns.begin(); it!=site_conns.end(); ++it) {
		if (!inVector<PseudoAtom>(*it, side1)) {
			side2.push_back(*it);
		}
	}

	// Verify that the pairings are self consistent
	if ((side1[0])->GetAngle(site, side1[1]) > 85) {
		obErrorLog.ThrowError(__FUNCTION__, "Trying to simplify a square-like four-connected pseudo atom", obWarning);
	}
	bool mismatched_sides = false;
	for (std::vector<PseudoAtom>::iterator it=site_conns.begin(); it!=site_conns.end(); ++it) {
		bool in_side1 = inVector<PseudoAtom>(*it, side1);
		if (in_side1) {
			if (!inVector<PseudoAtom>(minAngleNbor(site, *it), side1)) {
				mismatched_sides = true;
			}
		} else {  // If not in_side1, it should be part of side2
			if (!inVector<PseudoAtom>(minAngleNbor(site, *it), side2)) {
				mismatched_sides = true;
			}
		}
	}
	if (mismatched_sides) {
		obErrorLog.ThrowError(__FUNCTION__, "Inconsistent neighbors when assigning sides to the 4-c PA", obWarning);
	}

	// Start the PA split, now that sides (connections) have been assigned.
	VirtualMol side1_frag = VirtualMol(site);
	side1_frag.AddAtom(side1[0]);
	side1_frag.AddAtom(side1[1]);
	OBMol side1_mol = side1_frag.ToOBMol();
	vector3 side1_loc = getCentroid(&side1_mol, false);

	VirtualMol side2_frag = VirtualMol(site);
	side2_frag.AddAtom(side2[0]);
	side2_frag.AddAtom(side2[1]);
	OBMol side2_mol = side2_frag.ToOBMol();
	vector3 side2_loc = getCentroid(&side2_mol, false);

	// Use the original 4-c PA as one of the connections
	vector3 orig_site_loc = vector3(site->GetVector());  // copy the vector instead of using a reference
	site->SetVector(side1_loc);

	// And add a new 3-c site
	PseudoAtom new_3c = formAtom(&simplified_net, side2_loc, DEFAULT_ELEMENT);
	// no need to set act_to_pa since new_3c is empty
	pa_to_act[new_3c] = VirtualMol(orig_molp);  // empty: no atoms

	// Connect the new site to its two neighbors and the other 3-c
	vector3 side2_n0_loc = side2[0]->GetVector();
	vector3 side2_n1_loc = side2[1]->GetVector();
	ConnectAtoms(new_3c, site_x_to_nbors[side2[0]], &side2_n0_loc);
	ConnectAtoms(new_3c, site_x_to_nbors[side2[1]], &side2_n1_loc);
	ConnectAtoms(new_3c, site, &orig_site_loc);

	// Side 1 is no longer connected to side2
	DeleteConnection(side2[0]);
	DeleteConnection(side2[1]);

	return 1;
}


} // end namespace OpenBabel
