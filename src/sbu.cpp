// See https://openbabel.org/docs/dev/UseTheLibrary/CppExamples.html
// Get iterator help from http://openbabel.org/dev-api/group__main.shtml
// TODO: This works as a proof-of-concept, so then I'll need to fix the linker bonding,
// run tests, then analysis time!

#include <iostream>
#include <sstream>
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
void resetBonds(OBMol *mol);
bool subtractMols(OBMol *mol, OBMol *subtracted);
bool atomsEqual(const OBAtom &atom1, const OBAtom &atom2);


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
	obErrorLog.SetOutputLevel(obInfo);  // See also http://openbabel.org/wiki/Errors
	const bool display_full_smiles = false;
	const bool SAVE_NODES = true;
	bool DELETE_BONDS = true;
	char* filename;
	filename = argv[1];  // TODO: Check usage later

	OBMol mol;
	if (!readCIF(&mol, filename)) {
		printf("Error reading file: %s", filename);
		exit(1);
	}
	OBMol orig_mol = mol;  // Copy original definition to another variable for later use

	// Find linkers by deleting bonds to metals
	std::vector<OBMol> fragments;
	if (DELETE_BONDS) {
		obErrorLog.ThrowError(__FUNCTION__, "Bond deletion enabled", obDebug);
		// Same spirit as OBMol::DeleteHydrogens()
		std::vector<OBBond*> delbonds;
		FOR_ATOMS_OF_MOL(a, mol) {
			if (isMetal(&*a)) {
				FOR_BONDS_OF_ATOM(b, *a) {
					if (!inVector(&*b, delbonds)) {
						// avoid double deletion in the case of metal-metal bonds
						delbonds.push_back(&*b);
					}
				}
			}
		}
		mol.BeginModify();
		std::stringstream deletionMsg;
		for (std::vector<OBBond*>::iterator itb = delbonds.begin(); itb != delbonds.end(); ++itb) {
			deletionMsg << "Deleting bond with atomic numbers ";
			deletionMsg << (*itb)->GetBeginAtom()->GetAtomicNum() << " and ";
			deletionMsg << (*itb)->GetEndAtom()->GetAtomicNum();
			deletionMsg << " (OBBond * " << static_cast<void *>(&*itb) << ")\n";
			mol.DeleteBond(*itb);  // FIXME: this correctly deletes the bonds, but the atomic valencies are messed up like my original version.
			// IDEA: Maybe just need to rerun the ConnectTheDots type bond perception on the new fragments?
			// That approach most works, though we should also consider preserving the connection point instead of outright deleting the bond
		}
		obErrorLog.ThrowError(__FUNCTION__, deletionMsg.str(), obDebug);
		mol.EndModify();
	}


	fragments = mol.Separate();

	OBConversion obconv;
	obconv.SetOutFormat("can");
	obconv.AddOption("i");  // Ignore SMILES chirality for now

	// Get a list of unique SMILES code
	std::vector<std::string> unique_smiles;
	for (std::vector<OBMol>::iterator it = fragments.begin(); it != fragments.end(); ++it) {
		OBMol canon = *it;
		resetBonds(&canon);
		std::string mol_smiles = obconv.WriteString(&canon);
		if (!inVector<std::string>(mol_smiles, unique_smiles)) {
			unique_smiles.push_back(mol_smiles);
		}
	}

	printFragments(unique_smiles);  // Prints individual fragments from deleting bonds to metals

	// Now try extracting the nonmetal components
	// FIXME: currently a lot of copy-paste from the previous loop as proof-of-concept
	OBMol nodes = orig_mol;  // TODO: consider renaming, or having separate fragments vars for SMILES
	OBMol linkers = orig_mol;  // Can't do this by additions, because we need the UC data, etc.
	nodes.BeginModify();
	linkers.BeginModify();
	std::stringstream nonmetalMsg;
	for (std::vector<OBMol>::iterator it = fragments.begin(); it != fragments.end(); ++it) {
		std::string mol_smiles = obconv.WriteString(&*it);
		if (it->NumAtoms() == 1) {
			// std::cout << "Found an oxygen!\n";
			nonmetalMsg << "Found a solitary atom with atomic number " << it->GetFirstAtom()->GetAtomicNum() << std::endl;
			subtractMols(&linkers, &*it);
		//} else if () {
		//	std::cout << "Hydroxyl group" << std::endl;
		} else {
			nonmetalMsg << "Deleting linker " << mol_smiles;
			// linkers += *it;
			subtractMols(&nodes, &*it);
		}
		// Check for metal, single oxygens, etc
		// For starters, this could probably be done with strings, then use actual atom types later
		// Even better, try considering all single atoms and hydroxyl species as nodes.
		// If it's not a linker, loop over atoms of the *it molecule and delete them from nodes.
		// Then just extract nodes as what's left, as suggested in group meeting
	}
	obErrorLog.ThrowError(__FUNCTION__, nonmetalMsg.str(), obDebug);
	nodes.EndModify();
	linkers.EndModify();

	std::string whole_smiles = obconv.WriteString(&orig_mol);
	if (display_full_smiles) {
		printf("Original molecule: %s", whole_smiles.c_str());
		whole_smiles = obconv.WriteString(&mol);
		printf("New molecule: %s", whole_smiles.c_str());  // WriteString already outputs a newline
		whole_smiles = obconv.WriteString(&nodes);
		printf("Nodes: %s", whole_smiles.c_str());  // WriteString already outputs a newline
	}
	whole_smiles = obconv.WriteString(&nodes);
	printf("Nodes: %s", whole_smiles.c_str());  // WriteString already outputs a newline
	// ALSO LINKERS HERE

	if (SAVE_NODES) {
		OBConversion node_conv;
		node_conv.SetOutFormat("cif");  // mmcif has extra, incompatible fields
		node_conv.WriteFile(&nodes, "Test/nodes.cif");
		OBConversion linker_conv;
		linker_conv.SetOutFormat("cif");
		linker_conv.WriteFile(&linkers, "Test/linkers.cif");
	}

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

void resetBonds(OBMol *mol) {
	// Resets bond orders and bond detection for molecular fragments
	// FIXME: consider destroying all the bonds and starting from scratch
	// Also see https://sourceforge.net/p/openbabel/mailman/message/6229244/
	mol->ConnectTheDots();
	mol->PerceiveBondOrders();
}

bool subtractMols(OBMol *mol, OBMol *subtracted) {
	// Subtracts all atoms from a second molecule present in the first
	// Returns true if successful, false if failure (extra atoms, etc)
	// Modifies *mol if true, reverts back to original *mol if false
	// This code assumes that the molecule is well-defined (no multiply-defined atoms)

	std::vector<OBAtom*> deleted;  // TODO: consider implementing as a queue
	FOR_ATOMS_OF_MOL(a2, *subtracted) {
		bool matched = false;
		FOR_ATOMS_OF_MOL(a1, *mol) {
			if (atomsEqual(*a1, *a2)) {
				deleted.push_back(&*a1);
				matched = true;
				break;
			}
		}
		if (matched == false) {  // *subtracted has an extra atom
			return false;  // *mol not modified
		}
	}

	// It's okay to delete the atoms, since they should hold different (cloned) pointers
	mol->BeginModify();
	for (std::vector<OBAtom*>::iterator it = deleted.begin(); it != deleted.end(); ++it) {
		mol->DeleteAtom(*it);
	}
	mol->EndModify();
	return true;
}

bool atomsEqual(const OBAtom &atom1, const OBAtom &atom2) {
	// Loosely based on the data copied during OBAtom::Duplicate
	// Atoms are the same if these properties (excluding bond-related) are equivalent
	return atom1.GetAtomicNum() == atom2.GetAtomicNum() &&
			atom1.GetVector().IsApprox(atom2.GetVector(), 1.0e-6) &&
			atom1.GetIsotope() == atom2.GetIsotope();
}

