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
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/obiter.h>
#include <openbabel/babelconfig.h>
#include "config_sbu.h"


using namespace OpenBabel;  // See http://openbabel.org/dev-api/namespaceOpenBabel.shtml

// Function prototypes
bool readCIF(OBMol* molp, std::string filepath, bool bond_orders = true);
void writeCIF(OBMol* molp, std::string filepath, bool write_bonds = true);
void copyMOF(OBMol *src, OBMol *dest);
void writeSystre(OBMol* molp, std::string filepath, int element_x = 0, bool write_centers = true);
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
int collapseSBU(OBMol *mol, OBMol *fragment, int element = 118, int conn_element = 0);
int collapseTwoConn(OBMol* net, int ignore_element = 0);
int collapseXX(OBMol *net, int element_x);
int simplifyLX(OBMol *net, const std::vector<int> &linker_elements, int element_x);
vector3 getCentroid(OBMol *fragment, bool weighted);
std::vector<int> makeVector(int a, int b, int c);
vector3 unwrapFracNear(vector3 new_loc, vector3 ref_loc, OBUnitCell *uc);
vector3 unwrapCartNear(vector3 new_loc, vector3 ref_loc, OBUnitCell *uc);
OBBond* formBond(OBMol *mol, OBAtom *begin, OBAtom *end, int order = 1);
OBAtom* formAtom(OBMol *mol, vector3 loc, int element);

/* Define global parameters for MOF decomposition */
// Atom type for connection sites.  Assigned to Te (52) for now.  Set to zero to disable.
const int X_CONN = 52;
// Max. degrees between redundant L-X bonds.  Oxalic acid (extreme case) is < 85 degrees.
const int LX_ANGLE_TOL = 66;  // don't make this a round number without justification


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
		std::vector<int> used_elements() {
			std::vector<int> elements;
			for (std::map<std::string,int>::iterator it = _mapping.begin(); it != _mapping.end(); ++it) {
				elements.push_back(it->second);
			}
			return elements;
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

	OBMol orig_mol;
	// Massively improving performance by skipping kekulization of the full MOF
	if (!readCIF(&orig_mol, filename, false)) {
		printf("Error reading file: %s", filename);
		exit(1);
	}

	// Strip all of the original CIF labels, so they don't interfere with the automatically generated labels in the output
	FOR_ATOMS_OF_MOL(a, orig_mol) {
		if (a->HasData("_atom_site_label")) {
			a->DeleteData("_atom_site_label");
		}
	}

	/* Copy original definition to another variable for later use.
	 * Earlier, we performed all copies at once to reduce the performance bottleneck.
	 * In that case, copying internally called SSSR through EndModify(), so this statement was ~60% of the total code walltime.
	 * When these statements were nested together, the compiler was smart enough to amortize the operation.
	 * The biggest performance boost came from avoiding these perception behaviors altogether, which is the objective of copyMOF.
	 * As a result, certain structures, like ToBaCCo MOF bct_sym_10_mc_10__L_12.cif, no longer require 2+ minutes.
	 */
	OBMol mol, nodes, linkers, simplified_net;
	copyMOF(&orig_mol, &mol);
	copyMOF(&orig_mol, &nodes);
	copyMOF(&orig_mol, &linkers);  // Can't do this by additions, because we need the UC data, etc.
	copyMOF(&orig_mol, &simplified_net);

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
			collapseSBU(&simplified_net, &*it, linker_conv.key(mol_smiles), X_CONN);
		}
	}
	obErrorLog.ThrowError(__FUNCTION__, nonmetalMsg.str(), obDebug);
	nodes.EndModify();
	linkers.EndModify();
	simplified_net.EndModify();
	writeCIF(&simplified_net, "Test/test_partial.cif");

	// Simplify all the node SBUs into single points.
	ElementGen node_conv(true);
	simplified_net.BeginModify();
	std::vector<OBMol> sep_nodes = nodes.Separate();
	for (std::vector<OBMol>::iterator it = sep_nodes.begin(); it != sep_nodes.end(); ++it) {
		std::string node_smiles = getSMILES(*it, obconv);
		// Only calculating the connection points on the linkers
		collapseSBU(&simplified_net, &*it, node_conv.key(node_smiles));
	}
	simplified_net.EndModify();

	// For now, just print the SMILES of the nodes and linkers.
	// Eventually, this may become the basis for MOFFLES.
	printFragments(uniqueSMILES(nodes.Separate(), obconv));
	printFragments(uniqueSMILES(linkers.Separate(), obconv));

	// Write out the decomposed and simplified MOF
	writeCIF(&nodes, "Test/nodes.cif");
	writeCIF(&linkers, "Test/linkers.cif");
	writeCIF(&orig_mol, "Test/orig_mol.cif");
	writeCIF(&simplified_net, "Test/condensed_linkers.cif");
	writeFragmentKeys(node_conv.get_map(), linker_conv.get_map(), "Test/keys_for_condensed_linkers.txt");

	int simplifications;
	do {
		simplifications = 0;
		simplifications += collapseTwoConn(&simplified_net, X_CONN);
		simplifications += simplifyLX(&simplified_net, linker_conv.used_elements(), X_CONN);
		simplifications += collapseXX(&simplified_net, X_CONN);
	} while(simplifications);

	writeCIF(&simplified_net, "Test/removed_two_conn_for_topology.cif");
	writeSystre(&simplified_net, "Test/topology.cgd", X_CONN);

	return(0);
}


bool readCIF(OBMol* molp, std::string filepath, bool bond_orders) {
	// Read the first distinguished molecule from a CIF file
	// (TODO: check behavior of mmcif...)
	OBConversion obconversion;
	obconversion.SetInFormat("mmcif");
	obconversion.AddOption("p", OBConversion::INOPTIONS);
	if (!bond_orders) {
		obconversion.AddOption("s", OBConversion::INOPTIONS);
	}
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

void copyMOF(OBMol *src, OBMol *dest) {
	// Efficiently duplicates a MOF OBMol with single bonds.
	// Avoids implicitly calling expensive ring detection code on the full MOF by
	// temporarily disabling perception routines and filling in placeholder data.
	bool src_hybridization = src->HasHybridizationPerceived();
	if (!src_hybridization) {
		src->SetHybridizationPerceived();
		// Since bond orders aren't defined yet by OBMol::PerceiveBondOrders, let's set them to a known value so the values are initialized
		FOR_ATOMS_OF_MOL(a, *src) {
			a->SetHyb(0);
		}
	}
	bool src_atom_types = src->HasAtomTypesPerceived();
	if (!src_atom_types) {
		src->SetAtomTypesPerceived();
		FOR_ATOMS_OF_MOL(a, *src) {
			a->SetType("");
		}
	}
	bool src_charges = src->HasPartialChargesPerceived();
	if (!src_charges) {
		src->SetPartialChargesPerceived();
		FOR_ATOMS_OF_MOL(a, *src) {
			a->SetPartialCharge(0.0);
		}
	}

	(*dest) = (*src);  // OBMol copy

	// Restore the original flags
	if (!src_hybridization) {
		src->UnsetFlag(OB_HYBRID_MOL);
		dest->UnsetFlag(OB_HYBRID_MOL);
	}
	if (!src_atom_types) {
		src->UnsetFlag(OB_ATOMTYPES_MOL);
		dest->UnsetFlag(OB_ATOMTYPES_MOL);
	}
	if (!src_charges) {
		src->UnsetPartialChargesPerceived();
		dest->UnsetPartialChargesPerceived();
	}
}

void writeSystre(OBMol* pmol, std::string filepath, int element_x, bool write_centers) {
	// Write the simplified molecule to Systre for topological determination
	// element_x defines the atomic number of a special two-connected pseudo-atom.
	// If defined and nonzero, use X's as the graph edges.  Otherwise, the OBAtoms are directly linked through their OBBonds.
	// Can also print the (optional) edge_center field

	std::ofstream ofs;
	ofs.open(filepath.c_str());

	OBUnitCell* uc = pmol->GetPeriodicLattice();

	// Write header for the molecule
	std::string indent = "  ";
	ofs << "# CGD file generated by Open Babel " << BABEL_VERSION << ", see http://openbabel.sf.net" << std::endl;
	ofs << "CRYSTAL" << std::endl;
	ofs << indent << "NAME " << pmol->GetTitle() << std::endl;
	ofs << indent << "GROUP P1" << std::endl;  // TODO: use test cases to verify that Open Babel converts everything to P1
	ofs << indent << "CELL "
		<< uc->GetA() << " "
		<< uc->GetB() << " "
		<< uc->GetC() << " "
		<< uc->GetAlpha() << " "
		<< uc->GetBeta() << " "
		<< uc->GetGamma() << std::endl;

	int current_node = 0;
	FOR_ATOMS_OF_MOL(a, *pmol) {
		if (a->GetAtomicNum() != element_x) {
			++current_node;
			vector3 frac_coords = uc->CartesianToFractional(a->GetVector());
			ofs << indent << "NODE " << current_node
				<< " " << a->GetValence()  // coordination of the atom
				<< " " << frac_coords[0]
				<< " " << frac_coords[1]
				<< " " << frac_coords[2]
				<< std::endl;
		}
	}

	if (element_x) {
		std::vector<OBAtom*> visited_x;  // Have to keep track of visited connectors since we have tiny unit cells containing M-X-X-M'
		FOR_ATOMS_OF_MOL(x, *pmol) {
			if (x->GetAtomicNum() == element_x && !inVector<OBAtom*>(&*x, visited_x)) {
				visited_x.push_back(&*x);
				std::vector<vector3> vertices;  // Unwrap vertex coordinates near the connector X
				vector3 x_coords = uc->CartesianToFractional(x->GetVector());
				FOR_NBORS_OF_ATOM(n, *x) {
					vector3 raw_n_coords = uc->CartesianToFractional(n->GetVector());
					vector3 unwrap_n_coords = unwrapFracNear(raw_n_coords, x_coords, uc);
					if (n->GetAtomicNum() == element_x) {  // X-X present, so find the next connector.
						visited_x.push_back(&*n);
						// FIXME: This doesn't work due to the nature of unwrapping
						// x1 and x2 are both super close to M, so it's just bouncing back and forth on both sides.
						// Likely, we'll have to fix this by incorporating the two-connected linker.
						// Will it be generalizable if we transform X-L-X into X-X-X as an intermediate?
						FOR_NBORS_OF_ATOM(n2, *n) {
							if (&*n2 != &*x) {  // we don't want to return to x1
								if (n2->GetAtomicNum() == element_x) {
									obErrorLog.ThrowError(__FUNCTION__, "Unexpected X-X-X connection in simplified net.", obWarning);
								}
								vector3 raw_v2_coords = uc->CartesianToFractional(n2->GetVector());
								vector3 unwrap_v2_coords = unwrapFracNear(raw_v2_coords, unwrap_n_coords, uc);
								vertices.push_back(unwrap_v2_coords);
							}
						}
					} else {  // Standard vertex case, M-X-M
						vertices.push_back(unwrap_n_coords);
					}

				}
				if (vertices.size() != 2) {
					obErrorLog.ThrowError(__FUNCTION__, "Inconsistency: found connector X without two bonds", obWarning);
				}
				ofs << indent << "EDGE  "
					<< vertices[0][0] << " " << vertices[0][1]  << " " << vertices[0][2] << "   "
					<< vertices[1][0] << " " << vertices[1][1]  << " " << vertices[1][2] << std::endl;
			}
		}

		if (write_centers) {
			FOR_ATOMS_OF_MOL(x, *pmol) {
				if (x->GetAtomicNum() == element_x) {
					vector3 frac_coords = uc->CartesianToFractional(x->GetVector());
					ofs << indent << "# EDGE_CENTER  "
						<< frac_coords[0] << " "
						<< frac_coords[1] << " "
						<< frac_coords[2] << std::endl;
				}
			}
		}

	} else {
		std::stringstream edge_centers;
		FOR_BONDS_OF_MOL(b, *pmol) {
			std::vector<int> bond_dir = b->GetPeriodicDirection();
			vector3 f_add(bond_dir[0], bond_dir[1], bond_dir[2]);
			vector3 begin = uc->CartesianToFractional(b->GetBeginAtom()->GetVector());
			// For the second atom, we need to copy it to the correct unit cell with f_add
			vector3 end = uc->CartesianToFractional(b->GetEndAtom()->GetVector()) + f_add;
			ofs << indent << "EDGE  "
				<< begin[0] << " " << begin[1] << " " << begin[2] << "   "
				<< end[0] << " " << end[1] << " " << end[2] << std::endl;
			// Do edge center calculation, since we already have the coordinates handy
			edge_centers << "# EDGE_CENTER  "
				<< (begin[0] + end[0]) / 2.0 << " "
				<< (begin[1] + end[1]) / 2.0 << " "
				<< (begin[2] + end[2]) / 2.0 << std::endl;
		}

		if (write_centers)
			ofs << edge_centers.str();
	}

	ofs << "END" << std::endl;
	ofs.close();
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

int collapseSBU(OBMol *mol, OBMol *fragment, int element, int conn_element) {
	// Simplifies *mol by combining all atoms from *fragment into a single pseudo-atom
	// with atomic number element, maintaining connections to existing atoms.
	// If conn_element != 0, then add that element as a spacer at the connection
	// site (like "X" or "Q" atoms in top-down crystal generators).
	// Returns the number of external bonds maintained.

	std::vector<OBAtom*> orig_sbu;
	FOR_ATOMS_OF_MOL(a, *fragment) {
		OBAtom* orig_atom = atomInOtherMol(&*a, mol);
		if (!orig_atom) {
			obErrorLog.ThrowError(__FUNCTION__, "Tried to delete fragment not present in molecule", obError);
			return 0;
		}
		orig_sbu.push_back(orig_atom);
	}

	std::set<std::vector<OBAtom*> > connections;  // <External atom, internal atom bonded to it>
	for (std::vector<OBAtom*>::iterator it = orig_sbu.begin(); it != orig_sbu.end(); ++it) {
		FOR_NBORS_OF_ATOM(np, *it) {
			if (!inVector<OBAtom*>(&*np, orig_sbu)) {  // External atom: not in fragment
				std::vector<OBAtom*> new_connection;
				new_connection.push_back(&*np);
				new_connection.push_back(*it);
				connections.insert(new_connection);  // only push a connection once, as guaranteed by a set
			}
		}
	}

	vector3 centroid = getCentroid(fragment, false);

	mol->BeginModify();

	OBAtom* pseudo_atom = formAtom(mol, centroid, element);
	if (conn_element == 0) {  // Not using an "X" psuedo-atom
		std::vector<OBAtom*> new_external_conn;
		for (std::set<std::vector<OBAtom*> >::iterator it = connections.begin(); it != connections.end(); ++it) {
			OBAtom* external_atom = (*it)[0];
			if (!inVector<OBAtom*>(external_atom, new_external_conn)) {  // Only form external bonds once
				formBond(mol, pseudo_atom, external_atom, 1);
				new_external_conn.push_back(external_atom);
			}
		}
	} else {
		for (std::set<std::vector<OBAtom*> >::iterator it = connections.begin(); it != connections.end(); ++it) {
			// Put the connection point 1/3 of the way between the centroid and the connection point of the internal atom (e.g. COO).
			// In a simplified M-X-X-M' system, this will have the convenient property of being mostly equidistant.
			// Note: many top-down MOF generators instead place the connection point halfway on the node and linker bond.
			OBUnitCell* lattice = mol->GetPeriodicLattice();
			OBAtom* external_atom = (*it)[0];
			OBAtom* internal_atom = (*it)[1];

			vector3 internal_atom_loc = unwrapCartNear(internal_atom->GetVector(), centroid, lattice);
			vector3 conn_loc = lattice->WrapCartesianCoordinate((2.0*centroid + internal_atom_loc) / 3.0);

			OBAtom* conn_atom = formAtom(mol, conn_loc, conn_element);
			formBond(mol, conn_atom, pseudo_atom, 1);  // Connect to internal
			formBond(mol, conn_atom, external_atom, 1);
		}
	}

	// Delete SBU from the original molecule
	subtractMols(mol, fragment);

	mol->EndModify();

	return connections.size();
}

int collapseTwoConn(OBMol *net, int ignore_element) {
	// Collapses two-connected nodes into edges to simplify the topology
	// Returns the number of nodes deleted from the network

	int simplifications = 0;

	net->BeginModify();
	std::vector<OBAtom*> to_delete;
	FOR_ATOMS_OF_MOL(a, *net) {
		if (a->GetValence() == 2 && a->GetAtomicNum() != ignore_element) {
			std::vector<OBAtom*> nbors;
			FOR_NBORS_OF_ATOM(n, *a) {
				nbors.push_back(&*n);
			}

			formBond(net, nbors[0], nbors[1]);
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

int collapseXX(OBMol *net, int element_x) {
	// Simplify X-X bonds in the simplified net with their midpoint
	// FIXME: bond drawings in CIF are inconsistent due to ambiguity of unit cells for coordinates at 0.0
	// Returns the number of X-X bonds simplified
	int simplifications = 0;
	OBUnitCell* lattice = net->GetPeriodicLattice();

	net->BeginModify();
	bool new_simplification = true;
	while (new_simplification) {
		new_simplification = false;
		FOR_ATOMS_OF_MOL(a, *net) {
			OBAtom* x1 = &*a;
			if (x1->GetAtomicNum() == element_x) {  // Found the first X
				OBAtom* x2 = NULL;
				FOR_NBORS_OF_ATOM(n, *x1) {
					if (n->GetAtomicNum() == element_x) {
						x2 = &*n;
					}
				}
				if (x2) { // X-X exists, so simplify
					std::vector<OBAtom*> x_nbors;  // Get the neighbors of both X's
					FOR_NBORS_OF_ATOM(m, *x1) {
						if (&*m != x2) {
							x_nbors.push_back(&*m);
						}
					}
					FOR_NBORS_OF_ATOM(m, *x2) {
						if (&*m != x1) {
							x_nbors.push_back(&*m);
						}
					}
					if (x_nbors.size() != 2) {  // Each X should have one additional neighbor
						std::stringstream nborMsg;
						nborMsg << "y-X-X-y should be 4 atoms.  Found " << 2 + x_nbors.size() << std::endl;
						obErrorLog.ThrowError(__FUNCTION__, nborMsg.str(), obInfo);
					}
					if (x_nbors.size() == 2 && x_nbors[0] == x_nbors[1]) {
						// x1 and x2 are bonded to the same psuedo-atom (M-x1-x2-M').
						// They're not redundant, because M' is often the periodic image of M.
						// Simplifying x1-x2 will cause an inconsistency with two M-X bonds, so skip these x's
						continue;
					}

					// Make a new atom at the X-X midpoint
					vector3 unwrapped2 = unwrapCartNear(x2->GetVector(), x1->GetVector(), lattice);
					vector3 midpoint = (x1->GetVector() + unwrapped2) / 2.0;
					OBAtom* mid_atom = formAtom(net, midpoint, element_x);
					// Form bonds between the new midpoint atom and neighbors of X-X
					for (std::vector<OBAtom*>::iterator it = x_nbors.begin(); it != x_nbors.end(); ++it) {
						formBond(net, mid_atom, *it, 1);
					}
					// Delete old X-X atoms
					net->DeleteAtom(x1);
					net->DeleteAtom(x2);

					// Break out of the FOR_ATOMS loop so we can start with a fresh, unmodified loop
					new_simplification = true;
					break;
				}
			}
		}
		if (new_simplification) {
			simplifications += 1;
		}
	}
	net->EndModify();

	return simplifications;
}

int simplifyLX(OBMol *net, const std::vector<int> &linker_elements, int element_x) {
	// Remove redundant L-X bonds connecting the same linkers and nodes in the same direction
	// TODO: Consider averaging the positions (midpoint) instead of selecting one of the two at random.
	// Returns the number of L-X bonds deleted

	std::vector<OBAtom*> to_delete;

	FOR_ATOMS_OF_MOL(L, *net) {  // Iterating over linkers
		if (inVector<int>(L->GetAtomicNum(), linker_elements)) {  // if linker atom
			std::vector<OBAtom*> connectors;
			std::map<OBAtom*,OBAtom*> metals;  // <Atom X, M of L-X-M>
			FOR_NBORS_OF_ATOM(n, *L) {
				if (n->GetAtomicNum() == element_x) {  // connection site
					connectors.push_back(&*n);
					FOR_NBORS_OF_ATOM(M, *n) {  // So we can compare the node connections later
						if (!inVector<int>(M->GetAtomicNum(), linker_elements)) {
							metals[&*n] = &*M;
						}
					}
				}
			}

			// Iterate through the bonded L-X's to find redundant pairs
			for (std::vector<OBAtom*>::iterator x1 = connectors.begin(); x1 != connectors.end(); ++x1) {
				for (std::vector<OBAtom*>::iterator x2 = connectors.begin(); x2 != connectors.end(); ++x2) {
					// If the two X's are within an angle tolerance (and not already scheduled for deletion), delete the redundant copy.
					// By merit of the FOR loop over L, we've already established that L is the same.  Also check M with the metals map.
					if ( *x1 != *x2
						&& !inVector<OBAtom*>(*x1, to_delete)
						&& !inVector<OBAtom*>(*x2, to_delete)
						&& metals[*x1] == metals[*x2]
						&& (net->GetAngle(*x1, &*L, *x2) < LX_ANGLE_TOL) )
					{
						to_delete.push_back(*x2);
					}
				}
			}
		}
	}

	// Delete the redundant X's (which will also remove its X-L and X-M bonds)
	net->BeginModify();
	for (std::vector<OBAtom*>::iterator it = to_delete.begin(); it != to_delete.end(); ++it) {
		net->DeleteAtom(*it);
	}
	net->EndModify();

	return to_delete.size();
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

vector3 unwrapFracNear(vector3 new_loc, vector3 ref_loc, OBUnitCell *uc) {
	// Unwraps periodic, fractional coordinates (atom, etc.) at new_loc to be close to the reference
	// ref_loc, so you don't have to think about crossing box boundaries, etc.
	// i.e. unwrapNear(<0.9, 0.2, 0.2>, <0.3, 0.9, 0.2>) -> <-0.1, 1.2, 0.2>
	vector3 bond_dir = uc->PBCFractionalDifference(new_loc, ref_loc);
	return ref_loc + bond_dir;
}

vector3 unwrapCartNear(vector3 new_loc, vector3 ref_loc, OBUnitCell *uc) {
	// Unwraps periodic, Cartesian coordinates (atom, etc.) at new_loc to be close to the reference
	// ref_loc, so you don't have to think about crossing box boundaries, etc.
	// Similar ideas as unwrapFracNear
	vector3 bond_dir = uc->PBCCartesianDifference(new_loc, ref_loc);
	return ref_loc + bond_dir;
}

OBBond* formBond(OBMol *mol, OBAtom *begin, OBAtom *end, int order) {
	// Makes a bond between two atoms, complete with the proper accounting
	// TODO: decide how to handle cases where the bond already exists.
	// Overwrite existing bonds?  Make it, as long as it's a different periodic direction??
	OBBond* pseudo_link = NULL;
	if (!mol->GetBond(begin, end)) {
		pseudo_link = mol->NewBond();
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
	return pseudo_link;
}

OBAtom* formAtom(OBMol *mol, vector3 loc, int element) {
	// Makes a new atom with a specified location and atomic number
	OBAtom* atom = mol->NewAtom();
	atom->SetVector(loc);
	atom->SetAtomicNum(element);
	atom->SetType(etab.GetName(element));
	return atom;
}
