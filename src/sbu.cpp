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
	obconv.SetOutFormat("can");
	std::vector<std::string> unique_smiles;

	// Print out each fragment and get a list of unique SMILES code
	for (std::vector<OBMol>::iterator it = fragments.begin(); it != fragments.end(); ++it) {
		// Based on http://openbabel.org/dev-api/group__main.shtml
		printf("%d\t", it->NumAtoms());
		std::string mol_smiles = obconv.WriteString(&*it);
		printf("%s", mol_smiles.c_str());  // WriteString already outputs a newline
		bool new_flag = true;
		for (std::vector<std::string>::iterator i2 = unique_smiles.begin(); i2 != unique_smiles.end(); ++i2) {
			if (*i2 == mol_smiles) {
				new_flag = false;
				break;
			}
		}
		if (new_flag) {
			printf("(New!)\n");
			unique_smiles.push_back(mol_smiles);
		} else {
			printf("(Already defined)\n");
		}
	}
	printf("\n\nFound %d fragments:\n", unique_smiles.size());
	for (std::vector<std::string>::iterator i2 = unique_smiles.begin(); i2 != unique_smiles.end(); ++i2) {
		printf(" %s", i2->c_str());
	}

	return(0);
}


bool readCIF(OBMol* molp, std::string filepath) {
	// Read the first distinguished molecule from a CIF file
	// (TODO: check behavior of mmcif...)
	OBConversion obconversion;
	obconversion.SetInFormat("mmcif");
	return obconversion.ReadFile(molp, filepath);
}
