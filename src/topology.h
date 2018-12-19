/**********************************************************************
topology.h - Simplified topological net and related classes
***********************************************************************/

#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include <string>
#include <set>

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
	void SetRole(const std::string &possibly_new_role);
	bool DeleteRole(const std::string &role);
};


class Topology {
private:
	OBMol *orig_molp;
	OBMol simplified_net;
	VirtualMol connections;
	PseudoAtomMap atom_translator;  // map simplified PA to VirtualMol of orig atoms
	std::map<PseudoAtom, AtomRoles> role_translator;  // PA roles

	// Is a member of the simplified net a pseudo atom or connection?
	bool IsConnection(OBAtom* a);

public:
	//Topology() = delete;
	Topology(OBMol *parent_mol = NULL);
	OBMol* GetOrigMol() { return orig_molp; };
	VirtualMol GetOrigAtomsOfRole(const std::string &role);
	OBMol ToOBMol();
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
