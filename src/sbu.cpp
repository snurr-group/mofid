// See https://openbabel.org/docs/dev/UseTheLibrary/CppExamples.html

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>


/* #define filename "../Data/RingCIFs/DOTSOV_clean.cif" */

using namespace OpenBabel;  // See http://openbabel.org/dev-api/namespaceOpenBabel.shtml

bool readCIF(OBMol* molp, std::string filepath);

int main(int argc, char* argv[])
{
	char* filename;
	filename = argv[1];  // TODO: Check usage later

	OBMol mol;
	if (!readCIF(&mol, filename)) {
		printf("Error reading file: %s", filename);
		exit(1);
	}

	std::vector<OBMol> fragments;
	fragments = mol.Separate();

	OBConversion obconv;
	obconv.SetOutFormat("SMI");
	for (std::vector<OBMol>::iterator it = fragments.begin(); it != fragments.end(); ++it) {
		// Based on http://openbabel.org/dev-api/group__main.shtml
		printf("%d\t", it->NumAtoms());
		printf("%s\n", obconv.WriteString(&*it).c_str());
	}

	// Get a list of unique SMILES codes
	//std::vector<std::string> unique_smiles;
	// TODO: etc.....

	return(0);
}


bool readCIF(OBMol* molp, std::string filepath) {
	// Read the first distinguished molecule from a CIF file
	// (TODO: check behavior of mmcif...)
	OBConversion obconversion;
	obconversion.SetInFormat("mmcif");
	return obconversion.ReadFile(molp, filepath);
}
