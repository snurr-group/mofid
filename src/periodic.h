/**********************************************************************
periodic.h - Utilities to extend analysis on periodic OBMol's
***********************************************************************/

#ifndef PERIODIC_H
#define PERIODIC_H

#include <openbabel/babelconfig.h>
#include <map>

namespace OpenBabel
{
// forward declarations
class OBMol;
class OBAtom;
class OBBond;
class vector3;
class OBUnitCell;


class int3 {
// Simplified, int equivalent to vector3
public:
	int x;
	int y;
	int z;

	int3() {
		x=0; y=0; z=0;
	}
	int3(int a, int b, int c) {
		x=a; y=b; z=c;
	}
	// Get some of the operator functionality from vector3.cpp implementation
	bool operator== ( const int3 &other ) const {
		return ((x == other.x) && (y == other.y) && (z==other.z));
	}
	bool operator!= ( const int3 &other ) const {
		return !(*this == other);
	}
	/* Implementation from Open Babel, which does not allow assignment operations */
	/*
	int operator[] ( unsigned int i ) const {
		if (i == 0) { return x; }
		if (i == 1) { return y; }
		return z;  // otherwise
	}
	*/
	// In order to get the bracket assignment operator, we need to return a reference
	int& operator[] ( unsigned int i ) {
		if (i == 0) { return x; }
		if (i == 1) { return y; }
		return z;  // otherwise
	}
};


// Mapping of an atom to its relative unit cell/image, for an unwrapped molecule
typedef std::map<OBAtom*, int3> UCMap;

OBUnitCell* getPeriodicLattice(OBMol *mol);
bool isPeriodicChain(OBMol *mol);
int3 GetPeriodicDirection(OBBond *bond);
UCMap unwrapFragmentUC(OBMol *fragment, bool allow_rod = false, bool warn_rod = true);
bool unwrapFragmentMol(OBMol* fragment);
vector3 getCentroid(OBMol *fragment, bool weighted);
vector3 getMidpoint(OBAtom* a1, OBAtom* a2, bool weighted = false);
OBAtom* minAngleNbor(OBAtom* base, OBAtom* first_conn);

} // end namespace OpenBabel

#endif // PERIODIC_H

//! \file periodic.h
//! \brief periodic.h - Utilities to extend analysis on periodic OBMol's
