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
OBBond* formBond(OBMol *mol, OBAtom *begin, OBAtom *end, int order = 1);
OBAtom* formAtom(OBMol *mol, vector3 loc, int element);

} // end namespace OpenBabel

#endif // OB_DETAILS_H

//! \file obdetails.h
//! \brief Convenience functions to simplify use of Open Babel
