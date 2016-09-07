// See https://openbabel.org/docs/dev/UseTheLibrary/CppExamples.html
// Get iterator help from http://openbabel.org/dev-api/group__main.shtml

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>


using namespace OpenBabel;  // See http://openbabel.org/dev-api/namespaceOpenBabel.shtml

// Function prototypes
bool readCIF(OBMol* molp, std::string filepath);
void printFragments(const std::vector<std::string> &unique_smiles, bool display_number_of_fragments = false);
bool stringInVector(const std::string &element, const std::vector<std::string> &vec);


int main(int argc, char* argv[])
{
	// obErrorLog.SetOutputLevel(obDebug);  // See http://openbabel.org/wiki/Errors
	const bool display_full_smiles = false;
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
	obconv.AddOption("i");  // Ignore SMILES chirality for now

	if (display_full_smiles) {
		std::string whole_smiles = obconv.WriteString(&mol);
		printf("Whole molecule: %s", whole_smiles.c_str());  // WriteString already outputs a newline
	}

	// Get a list of unique SMILES code
	std::vector<std::string> unique_smiles;
	for (std::vector<OBMol>::iterator it = fragments.begin(); it != fragments.end(); ++it) {
		std::string mol_smiles = obconv.WriteString(&*it);
		if (!stringInVector(mol_smiles, unique_smiles)) {
			unique_smiles.push_back(mol_smiles);
		}
	}

	printFragments(unique_smiles);

	return(0);
}


bool readCIF(OBMol* molp, std::string filepath) {
	// Read the first distinguished molecule from a CIF file
	// (TODO: check behavior of mmcif...)
	OBConversion obconversion;
	obconversion.SetInFormat("mmcif");
	obconversion.AddOption("p", OBConversion::INOPTIONS);
	// Can disable bond detection as a diagnostic:
	// obconversion.AddOption("s", OBConversion::INOPTIONS);
	return obconversion.ReadFile(molp, filepath);
}

void printFragments(const std::vector<std::string> &unique_smiles, bool display_number_of_fragments) {
	// Print a list of fragments, optionally with a count
	if (display_number_of_fragments) {
		printf("\n\nFound %d fragments:\n", unique_smiles.size());
	}

	// Use a const_iterator since we're not modifying the vector: http://stackoverflow.com/questions/4890497/how-do-i-iterate-over-a-constant-vector
	for (std::vector<std::string>::const_iterator i2 = unique_smiles.begin(); i2 != unique_smiles.end(); ++i2) {
		printf("%s", i2->c_str());
	}
}

// template<class T>  // Consider using templates later, though there might be linker complications: https://isocpp.org/wiki/faq/templates#templates-defn-vs-decl
bool stringInVector(const std::string &element, const std::vector<std::string> &vec) {
	for (std::vector<std::string>::const_iterator it = vec.begin(); it != vec.end(); ++it) {
		if (*it == element) {
			return true;
		}
	}
	return false;
}
