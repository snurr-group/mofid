// See https://openbabel.org/docs/dev/UseTheLibrary/CppExamples.html
// Get iterator help from http://openbabel.org/dev-api/group__main.shtml

#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/obiter.h>


using namespace OpenBabel;  // See http://openbabel.org/dev-api/namespaceOpenBabel.shtml

// Function prototypes
bool readCIF(OBMol* molp, std::string filepath);
void printFragments(const std::vector<std::string> &unique_smiles, bool display_number_of_fragments = false);
bool isMetal(const OBAtom* atom);


template<typename T>  // WARNING: look out for linker complications: https://isocpp.org/wiki/faq/templates#templates-defn-vs-decl
bool inVector(const T &element, const std::vector<T> &vec) {
	// Test if an element is a given member of a vector
	if (std::find(vec.begin(), vec.end(), element) != vec.end()) {  // idiomatic C++
		return true;
	} else {
		return false;
	}
}


int main(int argc, char* argv[])
{
	// obErrorLog.SetOutputLevel(obDebug);  // See http://openbabel.org/wiki/Errors
	const bool display_full_smiles = false;
	bool DELETE_BONDS = true;
	char* filename;
	filename = argv[1];  // TODO: Check usage later

	OBMol mol;
	if (!readCIF(&mol, filename)) {
		printf("Error reading file: %s", filename);
		exit(1);
	}

	std::vector<OBMol> fragments;
	if (DELETE_BONDS) {
		printf("Bond deletion enabled\n");
		// Same spirit as OBMol::DeleteHydrogens()
		// std::set<OBBond*> delbonds;  // Can't actually use <set> and insert in STL because it's const
		std::vector<OBBond*> delbonds;
		FOR_ATOMS_OF_MOL(a, mol) {
			if (isMetal(&*a)) {
				FOR_BONDS_OF_ATOM(b, *a) {
					if (!inVector(&*b, delbonds)) {
						// avoid double deletion in the case of metal-metal bonds
						delbonds.push_back(&*b);
					}
					// See http://www.cplusplus.com/reference/set/set/insert/
				}
			}
		}
		mol.BeginModify();
		for (std::vector<OBBond*>::iterator itb = delbonds.begin(); itb != delbonds.end(); ++itb) {
			// TODO: Replace with OB debug statements throughout this module!
			printf("Deleting bond with atomic numbers %d and %d (OBBond* %p)\n",
					(*itb)->GetBeginAtom()->GetAtomicNum(),
					(*itb)->GetEndAtom()->GetAtomicNum(),
					static_cast<void *>(&*itb));
			mol.DeleteBond(*itb);  // FIXME: this correctly deletes the bonds, but the atomic valencies are messed up like my original version.
		}
		mol.EndModify();
	}


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
		if (!inVector<std::string>(mol_smiles, unique_smiles)) {
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

bool isMetal(const OBAtom* atom) {
	// Use the InChI definition of metals from https://jcheminf.springeropen.com/articles/10.1186/s13321-015-0068-4
	// Nonmetals are H, He, B, C, N, O, F, Ne, Si, P, S, Cl, Ar, Ge, As, Se, Br, Kr, Te, I, Xe, At, Rn
	const int NUM_NONMETALS = 23;
	unsigned int nonmetals[NUM_NONMETALS] = {1, 2, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18, 32, 33, 34, 35, 36, 52, 53, 54, 85, 86};

	unsigned int element = atom->GetAtomicNum();
	for (int i=0; i<NUM_NONMETALS; ++i) {
		if (nonmetals[i] == element) {
			return false;  // Atom is a nonmetal
		}
	}
	return true;  // Atom not classified as a nonmetal, so it's a metal
}
