/**********************************************************************
topology.h - Simplified topological net and related classes
***********************************************************************/

#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include <string>
#include <vector>
#include <set>
#include <tuple>

#include <openbabel/babelconfig.h>

#include "virtual_mol.h"
#include "pseudo_atom.h"

namespace OpenBabel
{
// forward declarations
class OBMol;
class OBAtom;

class AtomRoles {
// Roles of an orig_mol OBAtom in MOF simplification.
// Could potentially be re-implemented later as a bit vector.
// Example roles: SBU, organic, nodes, linkers, free_solvent, bound_solvent, solvent
private:
	std::set<std::string> _roles;
public:
	std::set<std::string> GetRoles() { return _roles; }
	bool HasRole(const std::string &test_role);
	void ClearRoles();
	void AddRole(const std::string &possibly_new_role);
	bool RemoveRole(const std::string &role);
};


class Connections {
private:
	OBMol *parent_net;
	std::map< PseudoAtom, std::pair<PseudoAtom, PseudoAtom> > conn2endpts;
	std::map< PseudoAtom, std::set<PseudoAtom> > endpt_nbors;
	std::map< PseudoAtom, std::set<PseudoAtom> > endpt_conns;
public:
	Connections(OBMol* parent = NULL);
	//bool CheckConsistency();
	void AddConn(PseudoAtom conn, PseudoAtom begin, PseudoAtom end);
	void RemoveConn(PseudoAtom conn);
	bool IsConn(PseudoAtom atom);
	AtomSet GetAtomConns(PseudoAtom atom);
	AtomSet GetAtomNeighbors(PseudoAtom atom);
	PseudoAtom GetConn(PseudoAtom begin, PseudoAtom end);
	AtomSet GetConnEndpoints(PseudoAtom conn);
	// Search for connection sites contained within a set of endpoints
	VirtualMol GetInternalConns(VirtualMol atoms);
};


class Topology {
private:
	const int DEFAULT_ELEMENT = 6;
	const int CONNECTION_ELEMENT = 118;  // Og for now

	OBMol *orig_molp;
	OBMol simplified_net;  // warning: see notes below about OBBonds
	Connections conns;
	VirtualMol deleted_atoms;
	PseudoAtomMap pa_to_act;  // map simplified PA to VirtualMol of orig atoms
	std::map<OBAtom*, AtomRoles> act_roles;  // roles of the original atoms
	std::map<OBAtom*, PseudoAtom> act_to_pa;  // where did the orig_mol atoms end up in the simplified net?

	// Is a member of the simplified net a pseudo atom or connection?
	bool IsConnection(PseudoAtom a);

public:
	//Topology() = delete;
	Topology(OBMol *parent_mol = NULL);
	OBMol* GetOrigMol() { return orig_molp; };
	VirtualMol GetOrigAtomsOfRole(const std::string &role);
	void SetRoleToAtoms(const std::string &role, VirtualMol atoms, bool val=true);
	int RemoveOrigAtoms(VirtualMol atoms);
	// Modify bonds using a custom connection-based routine rather than standard OBBonds:
	// Form a bond between two PseudoAtom's, taking care of all of the Connection accounting
	PseudoAtom ConnectAtoms(PseudoAtom begin, PseudoAtom end, vector3 *pos = NULL);
	void DeleteConnection(PseudoAtom conn);
	void DeleteConnection(PseudoAtom begin, PseudoAtom end);
	PseudoAtom CollapseOrigAtoms(VirtualMol atoms);
	OBMol ToOBMol();
	//something about deleting PA's?
};
// TODOs remaining
// How to handle PA connections?  Move/update SBU simplification methods from sbu.cpp?
// Same with topology simplification
// Export to Systre in systre.h or maybe an OBFormat
// Export topology as an OBMol with Og pseudoatoms?
// ToOBMol will require assigning pa_to_colors, based on unique SMILES etc.

} // end namespace OpenBabel
#endif // TOPOLOGY_H

//! \file topology.h
//! \brief topology.h - Simplified topological net and related classes
