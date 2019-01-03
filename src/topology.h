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


class ConnectionTable {
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
	VirtualMol deleted_atoms;
	PseudoAtomMap pa_to_act;  // map simplified PA to VirtualMol of orig atoms
	std::map<OBAtom*, std::string> pa_roles;  // roles of the simplified pseudoatoms
	std::map<OBAtom*, PseudoAtom> act_to_pa;  // where did the orig_mol atoms end up in the simplified net?
public:
	//Topology() = delete;
	Topology(OBMol *parent_mol = NULL);
	OBMol* GetOrigMol() { return orig_molp; };
	VirtualMol GetAtomsOfRole(const std::string &role);
	VirtualMol GetAtoms(bool include_conn=true);
	VirtualMol GetConnectors();
	bool AtomHasRole(PseudoAtom atom, const std::string &role);
	void SetRoleToAtom(const std::string &role, PseudoAtom atom);
	void SetRoleToAtoms(const std::string &role, VirtualMol atoms);
	std::string GetRoleFromAtom(PseudoAtom atom);
	int RemoveOrigAtoms(VirtualMol atoms);
	VirtualMol OrigToPseudo(VirtualMol orig_atoms);
	VirtualMol PseudoToOrig(VirtualMol pa_atoms);
	// Modify bonds using a custom connection-based routine rather than standard OBBonds:
	// Form a bond between two PseudoAtom's, taking care of all of the Connection accounting
	PseudoAtom ConnectAtoms(PseudoAtom begin, PseudoAtom end, vector3 *pos = NULL);
	void DeleteConnection(PseudoAtom conn);
	void DeleteAtomAndConns(PseudoAtom atom);
	PseudoAtom CollapseFragment(VirtualMol pa_fragment);
	void MergeAtomToAnother(PseudoAtom from, PseudoAtom to);
	OBMol FragmentToOBMolNoConn(VirtualMol pa_fragment);
	OBMol ToOBMol();
	void WriteSystre(const std::string &filepath, bool write_centers=true, bool simplify_two_conn=true);

	VirtualMol FragmentWithoutConns(VirtualMol fragment);
	VirtualMol FragmentWithIntConns(VirtualMol fragment);

	// Is a member of the simplified net a pseudo atom or connection?
	bool IsConnection(PseudoAtom a);
	//std::vector<VirtualMol> SeparatePeriodicRodIntoFragments(VirtualMol fragment);
	int SimplifyAxB();
	int SplitFourVertexIntoTwoThree(PseudoAtom site);
};


} // end namespace OpenBabel
#endif // TOPOLOGY_H

//! \file topology.h
//! \brief topology.h - Simplified topological net and related classes
