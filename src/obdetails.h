/**********************************************************************
obdetails.h - Convenience functions to simplify use of Open Babel
***********************************************************************/

#ifndef OB_DETAILS_H
#define OB_DETAILS_H

#include <openbabel/babelconfig.h>
#include <map>

namespace OpenBabel
{
// forward declarations
class OBAtom;
class OBBond;
class OBMol;
class vector3;

bool isMetal(const OBAtom* atom);
OBBond* formBond(OBMol *mol, OBAtom *begin, OBAtom *end, int order = 1);
OBAtom* formAtom(OBMol *mol, vector3 loc, int element);
int deleteBonds(OBMol *mol, bool only_metals = false);
bool subtractMols(OBMol *mol, OBMol *subtracted);
bool atomsEqual(const OBAtom &atom1, const OBAtom &atom2);
OBAtom* atomInOtherMol(OBAtom *atom, OBMol *mol);
bool isSubMol(OBMol *sub, OBMol *super);
std::map<int,int> getNumericFormula(OBMol *mol);

} // end namespace OpenBabel

#endif // OB_DETAILS_H

//! \file obdetails.h
//! \brief Convenience functions to simplify use of Open Babel
