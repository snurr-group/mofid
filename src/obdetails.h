/**********************************************************************
obdetails.h - Convenience functions to simplify use of Open Babel
***********************************************************************/

#ifndef OB_DETAILS_H
#define OB_DETAILS_H

#include <openbabel/babelconfig.h>

namespace OpenBabel
{
// forward declarations
class OBAtom;
class OBBond;
class OBMol;
class vector3;

// Class OBMol
bool isMetal(const OBAtom* atom);
OBBond* formBond(OBMol *mol, OBAtom *begin, OBAtom *end, int order = 1);
OBAtom* formAtom(OBMol *mol, vector3 loc, int element);
int deleteBonds(OBMol *mol, bool only_metals = false);

} // end namespace OpenBabel

#endif // OB_DETAILS_H

//! \file obdetails.h
//! \brief Convenience functions to simplify use of Open Babel
