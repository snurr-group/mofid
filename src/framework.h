/**********************************************************************
framework.h - Wrapper for OBMol with additional MOF/linker analyses
***********************************************************************/

#ifndef FRAMEWORK_H
#define FRAMEWORK_H

#include <openbabel/babelconfig.h>
#include <openbabel/generic.h>

#include <string>

namespace OpenBabel
{
// forward declarations
class OBMol;

extern bool COPY_ALL_CIFS_TO_PDB;

struct MinimalAtom {
// Contains the critial information about OBAtom's, to avoid accidentally copying perceived properties, etc.
	vector3 loc;
	int element;
	bool is_paddlewheel;
};

bool importCIF(OBMol* molp, std::string filepath, bool bond_orders = true, bool makeP1 = true);
void writeCIF(OBMol* molp, std::string filepath, bool write_bonds = true);
void writePDB(OBMol* orig_molp, std::string pdb_filepath, bool write_bonds = true, bool exclude_pbc = true);
OBMol initMOFwithUC(OBMol *orig_in_uc);
void copyMOF(OBMol *src, OBMol *dest);

void resetBonds(OBMol *mol);
void detectSingleBonds(OBMol *mol, double skin = 0.45, bool only_override_oxygen = true);
bool normalizeCharges(OBMol *mol);
bool detectPaddlewheels(OBMol *mol);

} // end namespace OpenBabel
#endif // FRAMEWORK_H

//! \file framework.h
//! \brief framework.h - Wrapper for OBMol with additional MOF/linker analyses
