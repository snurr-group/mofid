/**********************************************************************
virtual_mol.h - Collection of OBAtom*, plus adapters for OBMol, etc.
***********************************************************************/

#ifndef VIRTUAL_MOL_H
#define VIRTUAL_MOL_H

#include <string>
#include <set>

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>

namespace OpenBabel
{

class VirtualMol;

class VirtualMol {
private:
	std::set<OBAtom*> _atoms;
	OBMol *_parent_mol;
public:
	// VirtualMol() = delete;  // does not make sense when there's a default parameter below.
	// Note: the map STL container requires a no-argument constructor in case map[key] has an unknown key.
	// For more info, see https://stackoverflow.com/questions/695645/why-does-the-c-map-type-argument-require-an-empty-constructor-when-using
	VirtualMol(OBMol *parent = NULL);
	OBMol* GetParent();
	std::set<OBAtom*> GetAtoms();
	bool HasAtom(OBAtom *a);
	bool AddAtom(OBAtom *a);
	bool RemoveAtom(OBAtom *a);
	bool AddVirtualMol(VirtualMol addition);  // consider writing as operator+= or +
	OBMol ToOBMol(bool export_bonds = true, bool copy_bonds = true);
	// TODO: consider implementing SMILES in a parent class due to OBConv
	// std::string ToSmiles();}
};

} // end namespace OpenBabel
#endif // VIRTUAL_MOL_H

//! \file virtual_mol.h
//! \brief virtual_mol.h - Collection of OBAtom*, plus adapters for OBMol, etc.
