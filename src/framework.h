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

struct MinimalAtom {
	vector3 loc;
	int element;
	bool is_paddlewheel;
};

bool readCIF(OBMol* molp, std::string filepath, bool bond_orders = true, bool makeP1 = true);
void writeCIF(OBMol* molp, std::string filepath, bool write_bonds = true);
OBMol initMOF(OBMol *orig_in_uc);
void copyMOF(OBMol *src, OBMol *dest);

void resetBonds(OBMol *mol);
void detectSingleBonds(OBMol *mol, double skin = 0.45, bool only_override_oxygen = true);
bool normalizeCharges(OBMol *mol);
bool detectPaddlewheels(OBMol *mol);

} // end namespace OpenBabel
#endif // FRAMEWORK_H

//! \file framework.h
//! \brief framework.h - Wrapper for OBMol with additional MOF/linker analyses
