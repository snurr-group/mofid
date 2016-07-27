// See https://openbabel.org/docs/dev/UseTheLibrary/CppExamples.html
// Also see someone else's problem with compilation: http://forums.openbabel.org/Problem-for-compiling-the-example-files-td4656837.html
// In the current state, go to the build directory.  Use `cmake ..` then `make` to build the executable

#include <iostream>
#include <openbabel/obconversion.h>
using namespace std;

int main(int argc,char **argv)
{
  if(argc<3)
  {
    cout << "Usage: ProgrameName InputFileName OutputFileName\n";
    return 1;
  }

  ifstream ifs(argv[1]);
  if(!ifs)
  {
    cout << "Cannot open input file\n";
    return 1;
  }
  ofstream ofs(argv[2]);
  if(!ofs)
  {
    cout << "Cannot open output file\n";
    return 1;
  }
  OpenBabel::OBConversion conv(&ifs, &ofs);
  if(!conv.SetInAndOutFormats("CML","MOL"))
  {
    cout << "Formats not available\n";
    return 1;
  }
  int n = conv.Convert();
  cout << n << " molecules converted\n";

  return 0;
}
