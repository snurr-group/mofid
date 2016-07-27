// See https://openbabel.org/docs/dev/UseTheLibrary/CppExamples.html

#include <iostream>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>

#define FILENAME "../../Data/RingCIFs/DOTSOV_clean.cif"
#include <stdio.h>
#include <stdlib.h>
#include <openbabel/tokenst.h>

using namespace OpenBabel;

int main(int argc,char **argv)
{
  OBConversion obconversion;
  obconversion.SetInFormat("mmcif");
  OBMol mol;

  // From example at http://www.cplusplus.com/reference/cstdlib/getenv/
  char* pPath;
  pPath = getenv ("BABEL_DATADIR");
  if (pPath!=NULL) {
    printf ("The current BABEL_DATADIR is: %s\n",pPath);
    printf ("The contents of the directory are:\n");
    std::string cmd = "ls ";
    cmd += pPath;
    system(cmd.c_str());
    printf ("The directory separator character is %s\n", FILE_SEP_CHAR);
    // This finally works.  Just make sure that the file endings in the data are LF
    // Not sure if it's a bug in include/openbabel/lineend.h:169 or the symmetry code not calling it
    // From the looks of it, it seems like only molecule reads are done with that
    // method, not the data files.
  }
  
  bool notatend = obconversion.ReadFile(&mol, FILENAME);
  while (notatend)
  {
    std::cout << "Molecular Weight: " << mol.GetMolWt() << std::endl;

    mol.Clear();
    notatend = obconversion.Read(&mol);
  }

  return(0);
}
