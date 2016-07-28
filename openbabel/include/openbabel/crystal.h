/**********************************************************************
crystal.h - Handle crystals. Declarations of OBCryst.
        (the main header for Open Babel)

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

#ifndef OB_CRYST_H
#define OB_CRYST_H

#include <openbabel/babelconfig.h>

#include <openbabel/mol.h>

#include <math.h>
#include <float.h>

#include <vector>
#include <string>
#include <map>

namespace OpenBabel
{

  class OBCyrst;

  // class introduction in crystal.cpp
 class OBAPI OBCryst: public OBMol
  {
  protected:
    std::string                   _test;     	//!< Test for Ben
    OBUnitCell *                  _pCell;

  public:

    //! \name Initialization and data (re)size methods
    //@{
    //! Constructor
    OBCryst();
    //! Copy constructor, copies atoms,bonds and OBGenericData
    //OBCryst(const OBCryst &);
    //! Destructor
    //virtual ~OBCryst();
    //! Assignment, copies atoms,bonds and OBGenericData
    //OBMol &operator=(const OBMol &mol);
    //! Copies atoms and bonds but not OBGenericData
    //OBMol &operator+=(const OBMol &mol);

    //! Create a new OBBond in this molecule and ensure connections
    //! (e.g. OBBond::GetParent(). A new unique id will be assigned
    //! to this bond.
    //OBBond    *NewBond();
    //! Create a new OBBond in this molecule and ensure connections
    //! (e.g. OBBond::GetParent(). The @p id will be assigned to this
    //! bond.
    //OBBond    *NewBond(unsigned long id);
    //! Create a new OBResidue in this molecule and ensure connections.

    //! \name Data retrieval methods
    //@{
    //! \return the entire set of flags. (Internal use, mainly.)
    std::string  TestFunc(int flag)    { return((flag == 1) ? _test : "Not == 1"); }
    //@}

    //! \name Molecule utilities and perception methods
    //@{
    //! Adds single bonds based on atom proximity
    virtual void ConnectTheDots();
    //@}

  };

} // end namespace OpenBabel

#endif // OB_CRYST_H

//! \file crystal.h
//! \brief Handle crystals. Declarations of OBCryst AND OTHERS?  FIXME.
