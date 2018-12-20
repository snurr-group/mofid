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

const int CONNECTION_ELEMENT = 118;  // Og for now

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
	PseudoAtom FormConn(PseudoAtom conn, PseudoAtom begin, PseudoAtom end);
	PseudoAtom RemoveConn(PseudoAtom conn);
	bool IsConn(PseudoAtom atom);
	AtomSet GetAtomConns(PseudoAtom atom);
	AtomSet GetAtomNeighbors(PseudoAtom atom);
	AtomSet GetConnEndpoints(PseudoAtom conn);
	// Search for connection sites contained within a set of endpoints
	VirtualMol GetInternalConns(VirtualMol atoms);
};


class Topology {
private:
	OBMol *orig_molp;
	OBMol simplified_net;
	VirtualMol connections;
	VirtualMol deleted_atoms;
	PseudoAtomMap pa_to_act;  // map simplified PA to VirtualMol of orig atoms
	std::map<OBAtom*, AtomRoles> act_roles;  // roles of the original atoms
	std::map<OBAtom*, PseudoAtom> act_to_pa;  // where did the orig_mol atoms end up in the simplified net?

	// Is a member of the simplified net a pseudo atom or connection?
	bool IsConnection(OBAtom* a);
	std::vector<OBAtom*> GetConnectionPath(PseudoAtom conn);


public:
	//Topology() = delete;
	Topology(OBMol *parent_mol = NULL);
	OBMol* GetOrigMol() { return orig_molp; };
	VirtualMol GetOrigAtomsOfRole(const std::string &role);
	void SetRoleToAtoms(const std::string &role, VirtualMol atoms, bool val=true);
	int RemoveOrigAtoms(VirtualMol atoms);
	PseudoAtom CollapseOrigAtoms(VirtualMol atoms);
	OBMol ToOBMol();
	// todo: use FormConnection instead of an explicit FormBond
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
