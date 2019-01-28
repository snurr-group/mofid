/**********************************************************************
topology.h - Simplified topological net and related classes
***********************************************************************/

#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include <string>
#include <vector>
#include <set>
#include <utility>  // std::pair

#include <openbabel/babelconfig.h>

#include "virtual_mol.h"
#include "pseudo_atom.h"

namespace OpenBabel
{
// forward declarations
class OBMol;
class OBAtom;


const std::string DELETE_ORIG_ATOM_ERROR = "unexpected error";  // role for trying to delete PAs containing orig_mol atoms
const std::string ALL_DELETED_ORIG_ATOMS = "get all atoms";  // role for trying to delete PAs containing orig_mol atoms


class ConnectionTable {
// Handles accounting for connection pseudoatoms and their endpoints
private:
	OBMol *parent_net;
	std::map< PseudoAtom, std::pair<PseudoAtom, PseudoAtom> > conn2endpts;
	// Don't keep track of endpoint neighbors, since there may be multiple
	// connections between atoms 1 and 2 (e.g. different directions in IRMOF-1)
	// std::map< PseudoAtom, std::set<PseudoAtom> > endpt_nbors;
	std::map< PseudoAtom, std::set<PseudoAtom> > endpt_conns;
public:
	ConnectionTable(OBMol* parent = NULL);
	//bool CheckConsistency();
	void AddConn(PseudoAtom conn, PseudoAtom begin, PseudoAtom end);
	void RemoveConn(PseudoAtom conn);
	bool IsConn(PseudoAtom atom);
	AtomSet GetAtomConns(PseudoAtom endpt);
	bool HasNeighbor(PseudoAtom begin, PseudoAtom end);
	AtomSet GetConnEndpointSet(PseudoAtom conn);
	std::pair<PseudoAtom, PseudoAtom> GetConnEndpoints(PseudoAtom conn);
	PseudoAtom GetOtherEndpoint(PseudoAtom conn, PseudoAtom begin);
	// Search for connection sites contained within a set of endpoints
	VirtualMol GetInternalConns(VirtualMol atoms);
	// possibly consider a RemoveAtom(PseudoAtom endpt) to handle ConnectionTable accounting
	// but it wasn't necessary for Topology::DeleteAtom
};


class Topology {
// A simplified net, including explicit connections, of pseudoatoms which are initially copied
// from and mapped back to an original parent OBMol (the original MOF).

// Unlike an OBMol, PseudoAtoms are connected by an explicit connector OBAtom*, which allows for
// connecting the same two pseudoatoms in different directions (different UC's).

private:
	static const int DEFAULT_ELEMENT = 6;
	static const int CONNECTION_ELEMENT = 118;  // Og for now
	// Be careful with in-class constants.  They become trickier for non-integers:
	// http://www.stroustrup.com/bs_faq2.html#in-class

	OBMol *orig_molp;
	OBMol simplified_net;  // warning: see notes below about OBBonds
	// I'm wondering if simplified_net and related utilities should actually be a new class,
	// which would prevent inadvertently deleting bonds, etc.

	ConnectionTable conns;
	std::map<std::string, VirtualMol> deleted_atoms;  // atoms "deleted" from orig_molp in the simplified net
	PseudoAtomMap pa_to_act;  // map simplified PA to VirtualMol of orig atoms
	std::map<OBAtom*, std::string> pa_roles;  // roles of the simplified pseudoatoms
	std::map<OBAtom*, PseudoAtom> act_to_pa;  // where did the orig_mol atoms end up in the simplified net?

	// The complicated constructor makes a copy constructor nontrivial (and it's not currently being used).
	// Besides Wikipedia, here's another good overview: https://en.cppreference.com/w/cpp/language/rule_of_three
	Topology(const Topology& other);  // delete the copy constructor unless we need it and define it explicitly
	Topology& operator=(const Topology&);  // also copy assignment

public:
	//Topology() = delete;
	Topology(OBMol *parent_mol = NULL);
	OBMol* GetOrigMol() { return orig_molp; };

	// Manipulating/querying the roles of pseudoatoms
	bool AtomHasRole(PseudoAtom atom, const std::string &role);
	VirtualMol GetAtomsOfRole(const std::string &role);
	void SetRoleToAtom(const std::string &role, PseudoAtom atom);
	void SetRoleToAtoms(const std::string &role, VirtualMol atoms);
	std::string GetRoleFromAtom(PseudoAtom atom);

	// Exporting atoms and connections
	VirtualMol GetAtoms(bool include_conn=true);
	OBMol FragmentToOBMolNoConn(VirtualMol pa_fragment);
	VirtualMol GetDeletedOrigAtoms(const std::string &deletion_reason=ALL_DELETED_ORIG_ATOMS);
	VirtualMol GetConnectors();
	bool IsConnection(PseudoAtom a);
	PseudoAtom GetOtherEndpoint(PseudoAtom conn, PseudoAtom begin);

	// Conversions between the original and simplified nets
	VirtualMol OrigToPseudo(VirtualMol orig_atoms);
	VirtualMol PseudoToOrig(VirtualMol pa_atoms);
	// or removing/adding the explicit connection pseudoatoms
	VirtualMol FragmentWithoutConns(VirtualMol fragment);
	VirtualMol FragmentWithIntConns(VirtualMol fragment);

	// Manipulate connections (replaces direct manipulation of OBBond's)
	PseudoAtom ConnectAtoms(PseudoAtom begin, PseudoAtom end, vector3 *pos = NULL);
	void DeleteConnection(PseudoAtom conn);
	void DeleteAtomAndConns(PseudoAtom atom, const std::string &role_for_orig_atoms=DELETE_ORIG_ATOM_ERROR);

	// Network manipulation and simplification
	PseudoAtom CollapseFragment(VirtualMol pa_fragment);
	void MergeAtomToAnother(PseudoAtom from, PseudoAtom to);
	int SimplifyAxB();
	int SplitFourVertexIntoTwoThree(PseudoAtom site);

	// Export the Topology class to other formats
	OBMol ToOBMol();
	void ToSimplifiedCIF(const std::string &filename);
	void WriteSystre(const std::string &filepath, bool write_centers=true, bool simplify_two_conn=true);
};


} // end namespace OpenBabel
#endif // TOPOLOGY_H

//! \file topology.h
//! \brief topology.h - Simplified topological net and related classes
