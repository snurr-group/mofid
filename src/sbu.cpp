// See https://openbabel.org/docs/dev/UseTheLibrary/CppExamples.html
// Get iterator help from http://openbabel.org/dev-api/group__main.shtml
// Visualize with http://baoilleach.webfactional.com/site_media/blog/emscripten/openbabel/webdepict.html

// Instead of including code to display the full MOF SMILES, just use openbabel natively:
// obabel CIFFILE -ap -ocan

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <queue>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/obiter.h>
#include <openbabel/babelconfig.h>
#include "config_sbu.h"


using namespace OpenBabel;  // See http://openbabel.org/dev-api/namespaceOpenBabel.shtml

// Function prototypes
bool readCIF(OBMol* molp, std::string filepath);
void writeCIF(OBMol* molp, std::string filepath, bool write_bonds = true);
void writeSystre(OBMol* molp, std::string filepath);
void writeFragmentKeys(std::map<std::string,int> nodes, std::map<std::string,int> linkers, std::string filepath);
void printFragments(const std::vector<std::string> &unique_smiles);
std::string getSMILES(OBMol fragment, OBConversion obconv);
std::vector<std::string> uniqueSMILES(std::vector<OBMol> fragments, OBConversion obconv);
bool isMetal(const OBAtom* atom);
void resetBonds(OBMol *mol);
bool subtractMols(OBMol *mol, OBMol *subtracted);
int deleteBonds(OBMol *mol, bool only_metals = false);
bool atomsEqual(const OBAtom &atom1, const OBAtom &atom2);
OBAtom* atomInOtherMol(OBAtom *atom, OBMol *mol);
int collapseSBU(OBMol *mol, OBMol *fragment, int element = 118);
int collapseTwoConn(OBMol* net);
vector3 getCentroid(OBMol *fragment, bool weighted);
std::vector<int> makeVector(int a, int b, int c);
void formBond(OBMol *mol, OBAtom *begin, OBAtom *end, int order = 1);


class ElementGen
{
	protected:
		std::map<std::string,int> _mapping;
		std::queue<int> _elements;
		static const int _default_element = 118;  // Oganesson
		int _next() {
			int next = _elements.front();
			if (_elements.size() > 1) {
				_elements.pop();
			}
			return next;
		}
	public:
		ElementGen() {
			_elements.push(_default_element);
		}
		ElementGen(bool node) {
			if (node) {  // node defaults
				_elements.push(40);
				_elements.push(30);
				_elements.push(31);
				_elements.push(118);
				_elements.push(117);
			} else {  // linker defaults
				_elements.push(8);
				_elements.push(7);
				_elements.push(6);
				_elements.push(5);
			}
		}
		int key(std::string smiles) {
			if (_mapping.find(smiles) == _mapping.end()) {
				int new_key = _next();
				_mapping[smiles] = new_key;
				return new_key;
			}
			return _mapping[smiles];
		}
		std::map<std::string,int> get_map() {
			return _mapping;
		}
};


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
	char* filename;
	filename = argv[1];  // TODO: Check usage later

	// Set up the babel data directory to use a local copy customized for MOFs
	// (instead of system-wide Open Babel data)
	std::stringstream dataMsg;
	dataMsg << "Using local Open Babel data saved in " << LOCAL_OB_DATADIR << std::endl;
	obErrorLog.ThrowError(__FUNCTION__, dataMsg.str(), obAuditMsg);
	// Use setenv instead of putenv, per advice about string copies vs. pointers: http://stackoverflow.com/questions/5873029/questions-about-putenv-and-setenv/5876818#5876818
	// This is similar to the approach of cryos/avogadro:main.cpp:127
	// Per my objective, this only sets the environment within the scope of the sbu.exe program
	setenv("BABEL_DATADIR", LOCAL_OB_DATADIR, 1);

	OBMol mol;
	if (!readCIF(&mol, filename)) {
		printf("Error reading file: %s", filename);
		exit(1);
	}
	OBMol orig_mol = mol;  // Copy original definition to another variable for later use

	// Find linkers by deleting bonds to metals
	std::vector<OBMol> fragments;
	deleteBonds(&mol, true);
	fragments = mol.Separate();

	OBConversion obconv;
	obconv.SetOutFormat("can");
	obconv.AddOption("i");  // Ignore SMILES chirality for now

	// Get a list of unique SMILES code
	std::vector<std::string> unique_smiles = uniqueSMILES(fragments, obconv);
	std::stringstream fragmentMsg;
	fragmentMsg << "Unique fragments detected:" << std::endl;
	for (std::vector<std::string>::const_iterator i2 = unique_smiles.begin(); i2 != unique_smiles.end(); ++i2) {
		fragmentMsg << *i2;
	}
	obErrorLog.ThrowError(__FUNCTION__, fragmentMsg.str(), obDebug);

	// Classify nodes and linkers based on composition.
	// Consider all single atoms and hydroxyl species as node building materials.
	OBMol nodes = orig_mol;
	OBMol linkers = orig_mol;  // Can't do this by additions, because we need the UC data, etc.
	OBMol simplified_net = orig_mol;
	nodes.BeginModify();
	linkers.BeginModify();
	simplified_net.BeginModify();
	std::stringstream nonmetalMsg;
	ElementGen linker_conv(false);
	for (std::vector<OBMol>::iterator it = fragments.begin(); it != fragments.end(); ++it) {
		std::string mol_smiles = getSMILES(*it, obconv);
		// printf(mol_smiles.c_str());
		if (it->NumAtoms() == 1) {
			nonmetalMsg << "Found a solitary atom with atomic number " << it->GetFirstAtom()->GetAtomicNum() << std::endl;
			subtractMols(&linkers, &*it);
		} else if (mol_smiles == "[OH]\t\n") {
			// Tab is generated by the SMILES vs. name part of the string
			nonmetalMsg << "Found a hydroxyl group" << std::endl;
			subtractMols(&linkers, &*it);
		} else if (mol_smiles == "O\t\n") {
			nonmetalMsg << "Found a bound water molecule" << std::endl;
			subtractMols(&linkers, &*it);
		} else if (mol_smiles == "[O]O[O]\t\n") {
			nonmetalMsg << "Found a central oxygen with coordinated solvent" << std::endl;
			subtractMols(&linkers, &*it);
		} else {
			nonmetalMsg << "Deleting linker " << mol_smiles;
			subtractMols(&nodes, &*it);
			collapseSBU(&simplified_net, &*it, linker_conv.key(mol_smiles));
		}
	}
	obErrorLog.ThrowError(__FUNCTION__, nonmetalMsg.str(), obDebug);
	nodes.EndModify();
	linkers.EndModify();
	simplified_net.EndModify();

	// Simplify all the node SBUs into single points.
	ElementGen node_conv(true);
	simplified_net.BeginModify();
	std::vector<OBMol> sep_nodes = nodes.Separate();
	for (std::vector<OBMol>::iterator it = sep_nodes.begin(); it != sep_nodes.end(); ++it) {
		std::string node_smiles = getSMILES(*it, obconv);
		collapseSBU(&simplified_net, &*it, node_conv.key(node_smiles));
	}
	simplified_net.EndModify();

	// For now, just print the SMILES of the nodes and linkers.
	// Eventually, this may become the basis for MOFFLES.
	printFragments(uniqueSMILES(nodes.Separate(), obconv));
	printFragments(uniqueSMILES(linkers.Separate(), obconv));

	const bool EXPORT_NODES = true;
	if (EXPORT_NODES) {
		writeCIF(&nodes, "Test/nodes.cif");
		writeCIF(&linkers, "Test/linkers.cif");
		writeCIF(&orig_mol, "Test/orig_mol.cif");
		writeCIF(&simplified_net, "Test/condensed_linkers.cif");
		writeFragmentKeys(node_conv.get_map(), linker_conv.get_map(), "Test/keys_for_condensed_linkers.txt");
		collapseTwoConn(&simplified_net);
		writeSystre(&simplified_net, "Test/topology.cgd");
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

void writeCIF(OBMol* molp, std::string filepath, bool write_bonds) {
	// Write a molecule to file
	OBConversion conv;
	conv.SetOutFormat("cif");  // mmcif has extra, incompatible fields
	if (write_bonds) {
		conv.AddOption("g");
	}
	conv.WriteFile(molp, filepath);
}

void writeSystre(OBMol* molp, std::string filepath) {
	// Write the simplified molecule to Systre for topological determination
	OBConversion conv;
	conv.SetOutFormat("cgd");
	conv.AddOption("c");  // Optional flag to print edge centers
	conv.WriteFile(molp, filepath);
}

void writeFragmentKeys(std::map<std::string,int> nodes, std::map<std::string,int> linkers, std::string filepath) {
	// Save fragment identities for condensed_linkers.cif
	std::string equal_line = "================";
	std::ofstream out_file;
	out_file.open(filepath.c_str());
	out_file << "Fragment identities for condensed_linkers.cif" << std::endl << std::endl;

	out_file << "Nodes" << std::endl << equal_line << std::endl;
	for (std::map<std::string,int>::iterator it=nodes.begin(); it!=nodes.end(); ++it) {
		out_file << etab.GetSymbol(it->second) << ": " << it->first;
	}
	out_file << std::endl << std::endl;

	out_file << "Linkers" << std::endl << equal_line << std::endl;
	for (std::map<std::string,int>::iterator it=linkers.begin(); it!=linkers.end(); ++it) {
		out_file << etab.GetSymbol(it->second) << ": " << it->first;
	}
	out_file << std::endl;

	out_file.close();
}

void printFragments(const std::vector<std::string> &unique_smiles) {
	// Write a list of fragments
	// Use a const_iterator since we're not modifying the vector: http://stackoverflow.com/questions/4890497/how-do-i-iterate-over-a-constant-vector
	// TODO: consider stripping out extraneous tabs, etc, here or elsewhere in the code.
	for (std::vector<std::string>::const_iterator i2 = unique_smiles.begin(); i2 != unique_smiles.end(); ++i2) {
		printf("%s", i2->c_str());
	}
}

std::string getSMILES(OBMol fragment, OBConversion obconv) {
	// Prints SMILES based on OBConversion parameters
	OBMol canon = fragment;
	resetBonds(&canon);
	return obconv.WriteString(&canon);
}

std::vector<std::string> uniqueSMILES(std::vector<OBMol> fragments, OBConversion obconv) {
	// Extracts list of SMILES for unique fragments
	std::vector<std::string> unique_smiles;
	for (std::vector<OBMol>::iterator it = fragments.begin(); it != fragments.end(); ++it) {
		std::string mol_smiles = getSMILES(*it, obconv);
		if (!inVector<std::string>(mol_smiles, unique_smiles)) {
			unique_smiles.push_back(mol_smiles);
		}
	}
	return unique_smiles;
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

	// TODO: consider destroying all the bonds and starting from scratch
	// Also see https://sourceforge.net/p/openbabel/mailman/message/6229244/
	// For whatever reason, deleting all the bonds doesn't work out (disconnected fragments, etc.).  Need to diagnose why that's the case.
	// It fails for porphyrins.  I am guessing the cause is that Separate() does not copy over UC information, which is causing problems with PBC.
	// deleteBonds(mol, false);

	mol->BeginModify();
	FOR_ATOMS_OF_MOL(a, *mol) {
		a->SetFormalCharge(0);  // Not a specific reason for doing this, but it doesn't seem to make a difference.
		a->SetSpinMultiplicity(0);  // Reset radicals so that linker SMILES are consistent.
		a->SetHyb(0);  // Also reset hybridization in case that is causing problems
	}
	mol->ConnectTheDots();
	mol->PerceiveBondOrders();
	mol->EndModify();
}

bool subtractMols(OBMol *mol, OBMol *subtracted) {
	// Subtracts all atoms from a second molecule present in the first
	// Returns true if successful, false if failure (extra atoms, etc)
	// Modifies *mol if true, reverts back to original *mol if false
	// This code assumes that the molecule is well-defined (no multiply-defined atoms)

	std::vector<OBAtom*> deleted;  // Consider implementing as a queue
	FOR_ATOMS_OF_MOL(a2, *subtracted) {
		OBAtom* match = atomInOtherMol(&*a2, mol);
		if (!match) {  // *subtracted has an extra atom
			return false;  // *mol not modified
		}
		deleted.push_back(&*match);
	}

	// It's okay to delete the atoms, since they should hold different (cloned) pointers
	mol->BeginModify();
	for (std::vector<OBAtom*>::iterator it = deleted.begin(); it != deleted.end(); ++it) {
		mol->DeleteAtom(*it);
	}
	mol->EndModify();
	return true;
}

int deleteBonds(OBMol *mol, bool only_metals) {
	// Deletes bonds in a molecule, optionally only deleting bonds to metal atoms.
	// Returns the number of bond deletions.
	std::vector<OBMol> fragments;
	obErrorLog.ThrowError(__FUNCTION__, "Bond deletion enabled", obDebug);
	// Same spirit as OBMol::DeleteHydrogens()
	std::vector<OBBond*> delbonds;
	FOR_ATOMS_OF_MOL(a, *mol) {
		if (isMetal(&*a) || !only_metals) {
			FOR_BONDS_OF_ATOM(b, *a) {
				if (!inVector(&*b, delbonds)) {
					// avoid double deletion in the case of metal-metal bonds
					delbonds.push_back(&*b);
				}
			}
		}
	}
	mol->BeginModify();
	std::stringstream deletionMsg;
	for (std::vector<OBBond*>::iterator itb = delbonds.begin(); itb != delbonds.end(); ++itb) {
		deletionMsg << "Deleting bond with atomic numbers ";
		deletionMsg << (*itb)->GetBeginAtom()->GetAtomicNum() << " and ";
		deletionMsg << (*itb)->GetEndAtom()->GetAtomicNum();
		deletionMsg << " (OBBond * " << static_cast<void *>(&*itb) << ")\n";
		mol->DeleteBond(*itb);  // We should also consider preserving the connection point instead of outright deleting the bond
	}
	obErrorLog.ThrowError(__FUNCTION__, deletionMsg.str(), obDebug);
	mol->EndModify();
	return delbonds.size();
}

bool atomsEqual(const OBAtom &atom1, const OBAtom &atom2) {
	// Loosely based on the data copied during OBAtom::Duplicate
	// Atoms are the same if these properties (excluding bond-related) are equivalent
	return atom1.GetAtomicNum() == atom2.GetAtomicNum() &&
			atom1.GetVector().IsApprox(atom2.GetVector(), 1.0e-6) &&
			atom1.GetIsotope() == atom2.GetIsotope();
}

OBAtom* atomInOtherMol(OBAtom *atom, OBMol *mol) {
	// Is an atom within a molecule?  (as defined by atomsEqual conditions)
	// If so, return the pointer of the copycat.  Otherwise, return NULL
	// This has the nice benefit that the usage is also compatible with bool
	// ( if (atomInOtherMol(atomp, molp)) {}  will work, too.)
	FOR_ATOMS_OF_MOL(a, *mol) {
		if (atomsEqual(*a, *atom)) {
			return &*a;
		}
	}
	return NULL;
}

int collapseSBU(OBMol *mol, OBMol *fragment, int element) {
	// Simplifies *mol by combining all atoms from *fragment into a single pseudo-atom
	// with atomic number element, maintaining connections to existing atoms.
	// Returns the number of external bonds maintained.

	std::vector<OBAtom*> orig_sbu;
	FOR_ATOMS_OF_MOL(a, *fragment) {
		orig_sbu.push_back(atomInOtherMol(&*a, mol));
	}

	std::vector<OBAtom*> connections;
	for (std::vector<OBAtom*>::iterator it = orig_sbu.begin(); it != orig_sbu.end(); ++it) {
		FOR_NBORS_OF_ATOM(np, *it) {
			if ((!inVector<OBAtom*>(&*np, orig_sbu))  // External atom: not in fragment
					&& (!inVector<OBAtom*>(&*np, connections))) {  // only push once
				connections.push_back(&*np);
			}
		}
	}

	vector3 centroid = getCentroid(fragment, false);  // FIXME: check flags, usage, etc.

	// Delete SBU from the original molecule
	if (!subtractMols(mol, fragment)) {
		obErrorLog.ThrowError(__FUNCTION__, "Tried to delete fragment not present in molecule", obError);
		return 0;
	}

	mol->BeginModify();  // FIXME: should this be here or earlier?  Also consider the error handling

	OBAtom* pseudo_atom = mol->NewAtom();
	pseudo_atom->SetVector(centroid);
	pseudo_atom->SetAtomicNum(element);
	pseudo_atom->SetType(etab.GetName(element));

	for (std::vector<OBAtom*>::iterator it = connections.begin(); it != connections.end(); ++it) {
		formBond(mol, pseudo_atom, *it, 1);  // don't need to dereference iterator since it's a vector of pointers
	}

	mol->EndModify();

	return connections.size();
}

int collapseTwoConn(OBMol* net) {
	// Collapses two-connected nodes into edges to simplify the topology
	// Returns the number of nodes deleted from the network
	// FIXME: BUGGY LOGIC: doesn't have a problem with forming multiple single bonds between the same simplified net.

	int simplifications = 0;

	net->BeginModify();
	std::vector<OBAtom*> to_delete;
	FOR_ATOMS_OF_MOL(a, *net) {
		if (a->GetValence() == 2) {
			std::vector<OBAtom*> nbors;
			FOR_NBORS_OF_ATOM(n, *a) {
				nbors.push_back(&*n);
			}

			formBond(net, nbors[0], nbors[1]); // FIXME: probably needs a bond direction flag as well
			to_delete.push_back(&*a);
			++simplifications;
		}
	}
	for (std::vector<OBAtom*>::iterator it = to_delete.begin(); it != to_delete.end(); ++it) {
		net->DeleteAtom(*it);
	}
	net->EndModify();

	return simplifications;
}

vector3 getCentroid(OBMol *fragment, bool weighted) {
	// Calculate the centroid of a fragment, optionally weighted by the atomic mass
	// (which would give the center of mass).  Consider periodicity as needed, and
	// wrap coordinates inside of [0,1].
	// Returns a vector3 of the centroid's coordinates.

	vector3 center(0.0, 0.0, 0.0);
	double total_weight = 0;

	if (!fragment->IsPeriodic()) {
		FOR_ATOMS_OF_MOL(a, fragment) {
			double weight = 1.0;
			if (weighted) {
				weight = a->GetAtomicMass();
			}
			center += weight * a->GetVector();
			total_weight += weight;
		}
		return (center / total_weight);
	}

	// Otherwise, for periodic systems, we should build the fragment up atom-by-atom
	std::queue<OBAtom*> to_visit;
	std::map<OBAtom*, std::vector<int> > unit_cells;
	OBUnitCell* lattice = fragment->GetPeriodicLattice();

	// Start at whichever atom is (randomly?) saved first
	// Note: atom arrays begin with 1 in OpenBabel, while bond arrays begin with 0.
	OBAtom* start_atom = fragment->GetAtom(1);
	to_visit.push(start_atom);
	unit_cells[start_atom] = makeVector(0, 0, 0);  // original unit cell
	// Traverse the molecular graph until we run out of atoms
	while (!to_visit.empty()) {
		OBAtom* current = to_visit.front();
		to_visit.pop();
		FOR_NBORS_OF_ATOM(nbr, current) {
			if (unit_cells.find(&*nbr) == unit_cells.end()) {  // Unvisited atom
				// TODO: Consider implementing GetPeriodicDirection in atom.cpp as well
				OBBond* nbr_bond = fragment->GetBond(current, &*nbr);
				std::vector<int> uc = nbr_bond->GetPeriodicDirection();
				if (nbr_bond->GetBeginAtom() == &*nbr) {  // opposite bond direction as expected
					uc = makeVector(-1*uc[0], -1*uc[1], -1*uc[2]);
				}

				std::vector<int> current_uc = unit_cells[current];
				uc = makeVector(current_uc[0] + uc[0], current_uc[1] + uc[1], current_uc[2] + uc[2]);
				unit_cells[&*nbr] = uc;

				// Make sure to visit the neighbor (and its neighbors, etc.)
				// Each atom will only be traversed once, since we've already added it to unit_cells
				to_visit.push(&*nbr);
			}
		}
	}

	// We're assuming distinct, connected fragment SBUs, unless the system is nonperiodic.
	if (fragment->NumAtoms() != unit_cells.size()) {
		obErrorLog.ThrowError(__FUNCTION__, "getCentroid assumes one connected fragment for periodic systems.  Behavior is undefined.", obWarning);
		return vector3(0.0, 0.0, 0.0);
	}

	// Now that we've "unwrapped" the molecular fragment, calculate the centroid
	for (std::map<OBAtom*, std::vector<int> >::iterator it=unit_cells.begin(); it!=unit_cells.end(); ++it) {
		double weight = 1.0;
		if (weighted) {
			weight = it->first->GetAtomicMass();
		}
		std::vector<int> uc_shift = it->second;  // first: second = key: value of a map/dict.
		vector3 uc_shift_frac(uc_shift[0], uc_shift[1], uc_shift[2]);  // Convert ints to doubles
		vector3 coord_shift = lattice->FractionalToCartesian(uc_shift_frac);
		center += weight * (it->first->GetVector() + coord_shift);
		total_weight += weight;
	}
	center /= total_weight;
	center = lattice->WrapCartesianCoordinate(center);  // Keep result in the [0,1] unit cell
	return center;
}

std::vector<int> makeVector(int a, int b, int c) {
	// Shortcut for initializing a length-3 STL vector of ints.
	// TODO: Consider making a new int3 class similar to vector3.
	std::vector<int> v;
	v.push_back(a);
	v.push_back(b);
	v.push_back(c);
	return v;
}

void formBond(OBMol *mol, OBAtom *begin, OBAtom *end, int order) {
	// Makes a bond between two atoms, complete with the proper accounting
	// TODO: decide how to handle cases where the bond already exists.
	// Overwrite existing bonds?  Make it, as long as it's a different periodic direction??
	if (!mol->GetBond(begin, end)) {
		OBBond* pseudo_link = mol->NewBond();
		pseudo_link->SetBegin(begin);
		pseudo_link->SetEnd(end);
		pseudo_link->SetBondOrder(order);
		// Per OBBuilder::Connect in builder.cpp, we need to also update the atoms' bonding.
		// Otherwise, our OBAtomAtomIter will not operate properly (bonds will not propagate back to the atoms).
		begin->AddBond(pseudo_link);
		end->AddBond(pseudo_link);
	} else {
		obErrorLog.ThrowError(__FUNCTION__, "Did not generate multiply-defined bond between two atoms.", obWarning);
	}
}
