/********************************************************************** 
crystal.cpp - Handle crystals and periodic molecules.

Copyright (C) 2016 by Benjamin Bucior

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#include <openbabel/babelconfig.h>

#include <openbabel/crystal.h>
#include <openbabel/mol.h>
#include <openbabel/rotamer.h>
#include <openbabel/phmodel.h>
#include <openbabel/atomclass.h>
#include <openbabel/bondtyper.h>
#include <openbabel/builder.h>
#include <openbabel/math/matrix3x3.h>

#include <openbabel/stereo/tetrahedral.h>
#include <openbabel/stereo/cistrans.h>

#include <sstream>
#include <set>

using namespace std;

namespace OpenBabel
{

  // extern bool SwabInt;

  /** \class OBCryst <openbabel/periodic.h>
      \brief Crystal class for periodic molecules

      This class is derived from the OBMol class.
      It is designed to make minimal changes from OBMol to store properties
      relevant for crystals with minimal changes to OBMol.  The changes
      include bond detection across periodic boundaries (e.g., mmcif format)
      and bond distances based on periodic boundary conditions.

      An OBCryst class can be declared:
      \code
      OBCryst mol;
      \endcode

      And used exactly in place as OBMol.  The only difference comes if the
      unit cell is defined FIXME WITH EXAMPLE.

      See OBMol documentation for more generic methods and usage.
  */

  //
  // OBCryst member functions
  //

  // Check GetAngle, GetTorsion, SetTorsion, GetTorsion, FindAngles, FindTorsions
  // Also operator=, +=, Clear, BeginModify
  // And creation of the new fancy bonds


  OBCryst::OBCryst()
  {
    _pCell = new OBUnitCell;
  }

  OBCryst::~OBCryst()
  {
    // Clean up extra resources
    delete _pCell;
  }


  /*! This method adds single bonds between all atoms
    closer than their combined atomic covalent radii,
    then "cleans up" making sure bonded atoms are not
    closer than 0.4A and the atom does not exceed its valence.
    It implements blue-obelisk:rebondFrom3DCoordinates.

    Subtly overrides the ConnectTheDots methods from the base
    OBMol class.

  */
  void OBCryst::ConnectTheDots(void)
  {
    cout << "Overridden";
    if (Empty())
      return;
    if (_dimension != 3) return; // not useful on non-3D structures

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::ConnectTheDots", obAuditMsg);

    int j,k,max;
    double maxrad = 0;
    bool unset = false;
    OBAtom *atom,*nbr;
    vector<OBAtom*>::iterator i;
    vector<pair<OBAtom*,double> > zsortedAtoms;
    vector<double> rad;
    vector<int> zsorted;
    vector<int> bondCount; // existing bonds (e.g., from residues in PDB)

    double *c = new double [NumAtoms()*3];
    rad.resize(_natoms);

    for (j = 0, atom = BeginAtom(i) ; atom ; atom = NextAtom(i), ++j)
      {
        (atom->GetVector()).Get(&c[j*3]);
        pair<OBAtom*,double> entry(atom, atom->GetVector().z());
        zsortedAtoms.push_back(entry);
        bondCount.push_back(atom->GetValence());
      }
    sort(zsortedAtoms.begin(), zsortedAtoms.end(), SortAtomZ);

    max = zsortedAtoms.size();

    for ( j = 0 ; j < max ; j++ )
      {
        atom   = zsortedAtoms[j].first;
        rad[j] = etab.GetCovalentRad(atom->GetAtomicNum());
        maxrad = std::max(rad[j],maxrad);
        zsorted.push_back(atom->GetIdx()-1);
      }

    int idx1, idx2;
    double d2,cutoff,zd;
    for (j = 0 ; j < max ; ++j)
      {
    	double maxcutoff = SQUARE(rad[j]+maxrad+0.45);
        idx1 = zsorted[j];
        for (k = j + 1 ; k < max ; k++ )
          {
            idx2 = zsorted[k];

            // bonded if closer than elemental Rcov + tolerance
            cutoff = SQUARE(rad[j] + rad[k] + 0.45);

            zd  = SQUARE(c[idx1*3+2] - c[idx2*3+2]);
            // bigger than max cutoff, which is determined using largest radius,
            // not the radius of k (which might be small, ie H, and cause an early  termination)
            // since we sort by z, anything beyond k will also fail
            if (zd > maxcutoff )
              break;

            d2  = SQUARE(c[idx1*3]   - c[idx2*3]);
            if (d2 > cutoff)
              continue; // x's bigger than cutoff
            d2 += SQUARE(c[idx1*3+1] - c[idx2*3+1]);
            if (d2 > cutoff)
              continue; // x^2 + y^2 bigger than cutoff
            d2 += zd;

            if (d2 > cutoff)
              continue;
            if (d2 < 0.16) // 0.4 * 0.4 = 0.16
              continue;

            atom = GetAtom(idx1+1);
            nbr  = GetAtom(idx2+1);

            if (atom->IsConnected(nbr))
              continue;

            AddBond(idx1+1,idx2+1,1);
          }
      }

    // If between BeginModify and EndModify, coord pointers are NULL
    // setup molecule to handle current coordinates

    if (_c == NULL)
      {
        _c = c;
        for (atom = BeginAtom(i);atom;atom = NextAtom(i))
          atom->SetCoordPtr(&_c);
        _vconf.push_back(c);
        unset = true;
      }

    // Cleanup -- delete long bonds that exceed max valence
    OBBond *maxbond, *bond;
    double maxlength;
    vector<OBBond*>::iterator l, m;
    int valCount;
    bool changed;
    BeginModify(); //prevent needless re-perception in DeleteBond
    for (atom = BeginAtom(i);atom;atom = NextAtom(i))
      {
        while (atom->BOSum() > static_cast<unsigned int>(etab.GetMaxBonds(atom->GetAtomicNum()))
               || atom->SmallestBondAngle() < 45.0)
          {
            bond = atom->BeginBond(l);
            maxbond = bond;
            // Fix from Liu Zhiguo 2008-01-26
            // loop past any bonds
            // which existed before ConnectTheDots was called
            // (e.g., from PDB resdata.txt)
            valCount = 0;
            while (valCount < bondCount[atom->GetIdx() - 1]) {
              bond = atom->NextBond(l);
              // timvdm: 2008-03-05
              // NextBond only returns NULL if the iterator l == _bonds.end().
              // This was casuing problems as follows:
              // NextBond = 0x????????
              // NextBond = 0x????????
              // NextBond = 0x????????
              // NextBond = 0x????????
              // NextBond = NULL  <-- this NULL was not detected
              // NextBond = 0x????????
              if (!bond) // so we add an additional check
                break;
              maxbond = bond;
              valCount++;
            }
            if (!bond) // no new bonds added for this atom, just skip it
              break;

            // delete bonds between hydrogens when over max valence
            if (atom->IsHydrogen())
              {
                m = l;
                changed = false;
                for (;bond;bond = atom->NextBond(m))
                  {
                    if (bond->GetNbrAtom(atom)->IsHydrogen())
                      {
                        DeleteBond(bond);
                        changed = true;
                        break;
                      }
                  }
                if (changed)
                  {
                    // bond deleted, reevaluate BOSum
                    continue;
                  }
                else
                  {
                    // reset to first new bond
                    bond = maxbond;
                  }
              }

            maxlength = maxbond->GetLength();
            for (bond = atom->NextBond(l);bond;bond = atom->NextBond(l))
              {
                if (bond->GetLength() > maxlength)
                  {
                    maxbond = bond;
                    maxlength = bond->GetLength();
                  }
              }
            DeleteBond(maxbond); // delete the new bond with the longest length
          }
      }
    EndModify();
    if (unset)
      {
        _c = NULL;
        for (atom = BeginAtom(i);atom;atom = NextAtom(i))
          atom->ClearCoordPtr();
        _vconf.resize(_vconf.size()-1);
      }

    if (_c != NULL)
      delete [] c;
  }

  /*! Deleted */
  /*void OBMol::PerceiveBondOrders()
  {
    if (Empty())
      return;
    if (_dimension != 3) return; // not useful on non-3D structures

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::PerceiveBondOrders", obAuditMsg);

    OBAtom *atom, *b, *c;
    vector3 v1, v2;
    double angle;//, dist1, dist2;
    vector<OBAtom*>::iterator i;
    vector<OBBond*>::iterator j;//,k;

    //  BeginModify();

    // Pass 1: Assign estimated hybridization based on avg. angles
    for (atom = BeginAtom(i);atom;atom = NextAtom(i))
      {
        angle = atom->AverageBondAngle();  // FIXME: But atom.cpp:1013???

        //        cout << atom->GetAtomicNum() << " " << angle << endl;*/


} // end namespace OpenBabel

//! \file crystal.cpp
//! \brief Handle crystals. Implementation of OBCryst derived from OBMol.
