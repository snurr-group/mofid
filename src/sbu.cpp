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
#include <vector>
#include <map>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/obiter.h>
#include <openbabel/babelconfig.h>
#include <openbabel/phmodel.h>
#include <openbabel/elements.h>

#include "config_sbu.h"
#include "invector.h"
#include "obdetails.h"

using namespace OpenBabel;  // See http://openbabel.org/dev-api/namespaceOpenBabel.shtml


// Connections to an atom or fragment, in the form of a set where each element
// is the vector of OBAtom pointers: <External atom, internal atom bonded to it>
typedef std::set<std::vector<OBAtom*> > ConnExtToInt;

typedef std::map<std::string,std::set<OBAtom*> > MapOfAtomVecs;
typedef std::map<OBAtom*, std::vector<int> > UCMap;


// Function prototypes
std::string analyzeMOF(std::string filename);
extern "C" void analyzeMOFc(const char *cifdata, char *analysis, int buflen);
extern "C" int SmilesToSVG(const char* smiles, int options, void* mbuf, unsigned int buflen);
bool readCIF(OBMol* molp, std::string filepath, bool bond_orders = true, bool makeP1 = true);
void writeCIF(OBMol* molp, std::string filepath, bool write_bonds = true);
OBMol initMOF(OBMol *orig_in_uc);
void copyMOF(OBMol *src, OBMol *dest);
void writeSystre(OBMol* molp, std::string filepath, int element_x = 0, bool write_centers = true);
void writeFragmentKeys(std::map<std::string,int> nodes, std::map<std::string,int> linkers, std::map<std::string,int> removed, int connector, std::string filepath);
std::string writeFragments(const std::vector<std::string> &unique_smiles);
std::string getSMILES(OBMol fragment, OBConversion obconv);
std::vector<std::string> uniqueSMILES(std::vector<OBMol> fragments, OBConversion obconv);
void resetBonds(OBMol *mol);
void detectSingleBonds(OBMol *mol, double skin = 0.45, bool only_override_oxygen = true);
bool subtractMols(OBMol *mol, OBMol *subtracted);
bool atomsEqual(const OBAtom &atom1, const OBAtom &atom2);
OBAtom* atomInOtherMol(OBAtom *atom, OBMol *mol);
bool isSubMol(OBMol *sub, OBMol *super);
std::map<int,int> getNumericFormula(OBMol *mol);
ConnExtToInt getLinksToExt(OBMol *mol, OBMol *fragment);
std::vector<OBAtom*> uniqueExtAtoms(OBMol *mol, OBMol *fragment);
MapOfAtomVecs neighborsOverConn(OBAtom *loc, int skip_element);
OBAtom* collapseSBU(OBMol *mol, OBMol *fragment, int element = 118, int conn_element = 0);
int collapseTwoConn(OBMol* net, int ignore_element = 0);
int collapseXX(OBMol *net, int element_x);
int simplifyLX(OBMol *net, const std::vector<int> &linker_elements, int element_x);
std::vector<OBMol> removeOneConn(OBMol *net, std::map<OBAtom*, OBMol*> pseudo_to_real, std::vector<int> allowed_elements, int connector);
int fourToTwoThree(OBMol *net, int X_CONN);
OBAtom* minAngleNbor(OBAtom* base, OBAtom* first_conn);
UCMap unwrapFragmentUC(OBMol *fragment, bool allow_rod = false, bool warn_rod = true);
bool unwrapFragmentMol(OBMol* fragment);
vector3 getCentroid(OBMol *fragment, bool weighted);
vector3 getMidpoint(OBAtom* a1, OBAtom* a2, bool weighted = false);
bool isPeriodicChain(OBMol *mol);
int sepPeriodicChains(OBMol *nodes);
std::vector<int> makeVector(int a, int b, int c);
bool normalizeCharges(OBMol *mol);
bool detectPaddlewheels(OBMol *mol);
std::vector<int> GetPeriodicDirection(OBBond *bond);
OBUnitCell* getPeriodicLattice(OBMol *mol);

/* Define global parameters for MOF decomposition */
// Atom type for connection sites.  Assigned to Te (52) for now.  Set to zero to disable.
const int X_CONN = 52;


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
		int remove_key(std::string smiles) {
			if (_mapping.find(smiles) == _mapping.end()) {
				return 0;
			}
			int old_element = _mapping[smiles];
			_mapping.erase(smiles);
			return old_element;
		}
};

struct MinimalAtom {
	vector3 loc;
	int element;
	bool is_paddlewheel;
};


int main(int argc, char* argv[])
{
	obErrorLog.SetOutputLevel(obWarning);  // See also http://openbabel.org/wiki/Errors
	char* filename;
	filename = argv[1];  // TODO: Check usage later

	// Set up the babel data directory to use a local copy customized for MOFs
	// (instead of system-wide Open Babel data)
	std::stringstream dataMsg;
	dataMsg << "Using local Open Babel data saved in " << LOCAL_OB_DATADIR << std::endl;
	obErrorLog.ThrowError(__FUNCTION__, dataMsg.str(), obAuditMsg);
	// Use setenv instead of putenv, per advice about string copies vs. pointers: http://stackoverflow.com/questions/5873029/questions-about-putenv-and-setenv/5876818#5876818
	// This is similar to the approach of cryos/avogadro:main.cpp:127
	// Per my objective, this only sets the environment within the scope of the sbu.exe program.
	// But Windows defines a separate _putenv, etc: https://stackoverflow.com/questions/17258029/c-setenv-undefined-identifier-in-visual-studio
#ifdef _WIN32
	_putenv_s("BABEL_DATADIR", LOCAL_OB_DATADIR);
#else
#ifndef __INTELLISENSE__ // Ignore setenv error in vscode:
	setenv("BABEL_DATADIR", LOCAL_OB_DATADIR, 1);
#endif  // vscode error workaround
#endif

	std::string mof_results = analyzeMOF(std::string(filename));
	if (mof_results == "") {  // No MOFs found
		return(1);
	} else {
		std::cout << mof_results;
		return(0);
	}
}

std::string analyzeMOF(std::string filename) {
	// Extract components of the MOFid
	// Reports nodes/linkers, number of nets found, and writes CIFs to Test/ directory

	std::stringstream analysis;
	OBMol orig_mol;
	// Massively improving performance by skipping kekulization of the full MOF
	if (!readCIF(&orig_mol, filename, false)) {
		std::cerr << "Error reading file: %s" << filename << std::endl;
		return "";
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
	OBMol mol, nodes, linkers, simplified_net, mof_fsr, mof_asr;
	copyMOF(&orig_mol, &mol);
	copyMOF(&orig_mol, &nodes);
	copyMOF(&orig_mol, &linkers);  // Can't do this by additions, because we need the UC data, etc.
	copyMOF(&orig_mol, &simplified_net);
	copyMOF(&orig_mol, &mof_fsr);
	copyMOF(&orig_mol, &mof_asr);
	OBMol free_solvent = initMOF(&orig_mol);
	OBMol bound_solvent = initMOF(&orig_mol);

	writeCIF(&orig_mol, "Test/orig_mol.cif");

	std::ofstream file_info;
	file_info.open("Test/mol_name.txt", std::ios::out | std::ios::trunc);
	if (file_info.is_open()) {
		file_info << filename << std::endl;
		file_info.close();
	}


	// Find linkers by deleting bonds to metals
	std::vector<OBMol> fragments;
	deleteBonds(&mol, true);
	fragments = mol.Separate();

	OBConversion obconv;
	// Universal SMILES:
	//obconv.SetOutFormat("smi");
	//obconv.AddOption("U");
	// InChI or InChIKey, with same flags as Universal SMILES:
	//obconv.SetOutFormat("inchi");
	//obconv.SetOutFormat("inchikey");
	//obconv.AddOption("F");
	//obconv.AddOption("M");
	// Open Babel canonical SMILES:
	obconv.SetOutFormat("can");
	obconv.AddOption("i");  // Ignore SMILES chirality for now

	// Classify nodes and linkers based on composition.
	// Consider all single atoms and hydroxyl species as node building materials.
	nodes.BeginModify();
	linkers.BeginModify();
	simplified_net.BeginModify();
	free_solvent.BeginModify();
	mof_fsr.BeginModify();
	mof_asr.BeginModify();
	std::stringstream nonmetalMsg;
	ElementGen linker_conv(false);
	std::map<OBAtom*, OBMol*> pseudo_map;  // <pseudo atom location, fragment used to simplify it>
	for (std::vector<OBMol>::iterator it = fragments.begin(); it != fragments.end(); ++it) {
		std::string mol_smiles = getSMILES(*it, obconv);
		// If str comparisons are required, include a "\t\n" in the proposed smiles
		bool all_oxygens = true;  // Also allow hydroxyls, etc.
		FOR_ATOMS_OF_MOL(a, *it){
			if (a->GetAtomicNum() != 8 && a->GetAtomicNum() != 1) {
				all_oxygens = false;
			}
		}

		if (uniqueExtAtoms(&simplified_net, &*it).size() == 0) {
			// Assume free solvents are organic (or lone metals), so they'd be isolated without any external connections
			nonmetalMsg << "Deleting free solvent " << mol_smiles;
			free_solvent += *it;
			subtractMols(&linkers, &*it);
			subtractMols(&nodes, &*it);
			subtractMols(&simplified_net, &*it);
			subtractMols(&mof_fsr, &*it);
			subtractMols(&mof_asr, &*it);
		} else if (it->NumAtoms() == 1) {
			nonmetalMsg << "Found a solitary atom with atomic number " << it->GetFirstAtom()->GetAtomicNum() << std::endl;
			subtractMols(&linkers, &*it);
		} else if (all_oxygens) {
			nonmetalMsg << "Found an oxygen species " << mol_smiles;
			subtractMols(&linkers, &*it);
		} else {
			nonmetalMsg << "Deleting linker " << mol_smiles;
			subtractMols(&nodes, &*it);
			OBAtom* pseudo_atom = collapseSBU(&simplified_net, &*it, linker_conv.key(mol_smiles), X_CONN);
			pseudo_map[pseudo_atom] = &*it;
		}
	}
	obErrorLog.ThrowError(__FUNCTION__, nonmetalMsg.str(), obDebug);
	nodes.EndModify();
	linkers.EndModify();
	simplified_net.EndModify();
	free_solvent.EndModify();
	mof_fsr.EndModify();
	mof_asr.EndModify();
	writeCIF(&simplified_net, "Test/test_partial.cif");

	// Simplify all the node SBUs into single points.
	ElementGen node_conv(true);
	simplified_net.BeginModify();

	bool mil_type_mof = false;
	if (sepPeriodicChains(&nodes)) {
		mil_type_mof = true;
	}

	std::vector<OBMol> sep_nodes = nodes.Separate();
	for (std::vector<OBMol>::iterator it = sep_nodes.begin(); it != sep_nodes.end(); ++it) {
		std::string node_smiles = getSMILES(*it, obconv);
		OBAtom* pseudo_node;
		if (node_smiles == "[O]\t\n") {  // These oxygen atoms behave more like linkers
			pseudo_node = collapseSBU(&simplified_net, &*it, node_conv.key(node_smiles), X_CONN);
		} else {
			// Don't need extra connection atoms from the nodes
			pseudo_node = collapseSBU(&simplified_net, &*it, node_conv.key(node_smiles));
		}
		pseudo_map[pseudo_node] = &*it;
	}
	simplified_net.EndModify();


	// Handle one-connected species, notably bound solvents and metal-containing ligands.
	// Consider making this a do-while loop while the one-connected ends exist?
	simplified_net.BeginModify();
	nodes.BeginModify();
	linkers.BeginModify();
	std::vector<OBAtom*> pseudo_deletions;  // Save deletions for the end to avoid invalidating the iterator over simplified_net
	FOR_ATOMS_OF_MOL(a, simplified_net) {
		std::set<OBAtom*> a_nbors = neighborsOverConn(&*a, X_CONN)["nbors"];
		if (a_nbors.size() == 1) {
			int point_type = a->GetAtomicNum();

			if (inVector<int>(point_type, node_conv.used_elements())) {  // Metal-containing linker
				if (X_CONN) {  // Simplify by removing the connection point, simplifying the problem to the base case
					std::set<OBAtom*> a_conns = neighborsOverConn(&*a, X_CONN)["skipped"];
					for (std::set<OBAtom*>::iterator it=a_conns.begin(); it!=a_conns.end(); ++it) {
						pseudo_deletions.push_back(*it);  // Delete the relevant connector atoms at the end
					}
					formBond(&simplified_net, *(a_nbors.begin()), &*a);
				}

				// Reclassify the node atoms as part of the linker
				OBMol* fragment_mol = pseudo_map[&*a];
				linkers += *fragment_mol;
				subtractMols(&nodes, fragment_mol);

				// Use the original CIF to reconnect the metal-containing linker
				ConnExtToInt orig_links = getLinksToExt(&orig_mol, fragment_mol);
				for (ConnExtToInt::iterator it=orig_links.begin(); it!=orig_links.end(); ++it) {
					OBAtom* node_pos = atomInOtherMol((*it)[1], &linkers);
					OBAtom* linker_pos = atomInOtherMol((*it)[0], &linkers);
					formBond(&linkers, node_pos, linker_pos, 1);  // All bonds in the original CIF are single bonds by default
				}

				// At the end, delete the one-connected "nodes"
				pseudo_deletions.push_back(&*a);
				pseudo_map.erase(&*a);

			} else if (inVector<int>(point_type, linker_conv.used_elements())) {
				// Bound ligands, such as capping agents or bound solvents for ASR removal
				// This will be handled below after the nets are fully simplified
				continue;
			} else if (point_type == X_CONN) {
				obErrorLog.ThrowError(__FUNCTION__, "Found singly-connected connection atom in the simplified net.", obError);
			} else {
				std::string err_msg =
					"Unexpected atom type " +
					std::string(OBElements::GetName(point_type)) +
					"in the simplified net.";
				obErrorLog.ThrowError(__FUNCTION__, err_msg, obError);
			}
		}
	}
	for (std::vector<OBAtom*>::iterator it=pseudo_deletions.begin(); it!=pseudo_deletions.end(); ++it) {
		simplified_net.DeleteAtom(*it);
	}
	// Self-consistency bookkeeping: update pseudo atom map
	std::vector<OBMol> split_linkers = linkers.Separate();
	for (std::map<OBAtom*, OBMol*>::iterator it=pseudo_map.begin(); it!=pseudo_map.end(); ++it) {
		for (std::vector<OBMol>::iterator it2=split_linkers.begin(); it2!=split_linkers.end(); ++it2) {
			if (isSubMol(it->second, &*it2)) {
				it->second = &*it2;  // use the updated linker fragment
				break;
			}
		}
	}
	// Also update pseudo-atom element types for the linkers (ignoring nodes and linkers)
	FOR_ATOMS_OF_MOL(a, simplified_net) {
		int a_element = a->GetAtomicNum();
		if (a_element != X_CONN && !inVector<int>(a_element, node_conv.used_elements())) {
			std::string a_smiles = getSMILES(*(pseudo_map[&*a]), obconv);
			a->SetAtomicNum(linker_conv.key(a_smiles));  // Looks for an old key, or assigns a new one if it's an updated linker
		}
	}
	simplified_net.EndModify();
	nodes.EndModify();
	linkers.EndModify();


	writeCIF(&simplified_net, "Test/condensed_linkers.cif");


	// Catenation: check that all interpenetrated nets contain identical components.
	// If X_CONN tends to be overly inconsistent, we could remove it from the formula and
	// resimplify the coefficients for just the nodes and linkers.
	// Note: is this test really necessary if we check for topology later?
	std::vector<OBMol> net_components = simplified_net.Separate();
	std::string base_formula = "";
	if (!net_components.size()) {
		obErrorLog.ThrowError(__FUNCTION__, "No MOFs found in the simplified net.", obError);
	} else {
		base_formula = net_components[0].GetFormula();
	}
	for (std::vector<OBMol>::iterator it = net_components.begin(); it != net_components.end(); ++it) {
		std::string component_formula = it->GetFormula();
		if (component_formula != base_formula) {
			std::string err_msg =
				"Inconsistency in catenated nets.  Net with simplified formula " +
				component_formula + " does not match " + base_formula;
			obErrorLog.ThrowError(__FUNCTION__, err_msg, obWarning);
		}
	}

	bound_solvent.BeginModify();
	mof_asr.BeginModify();
	linkers.BeginModify();
	int simplifications;
	do {
		simplifications = 0;
		simplifications += simplifyLX(&simplified_net, linker_conv.used_elements(), X_CONN);
		simplifications += collapseXX(&simplified_net, X_CONN);

		std::vector<OBMol> one_conn = removeOneConn(&simplified_net, pseudo_map, linker_conv.used_elements(), X_CONN);
		for (std::vector<OBMol>::iterator it=one_conn.begin(); it!=one_conn.end(); ++it) {
			bound_solvent += *it;
			subtractMols(&linkers, &*it);
		}
		simplifications += one_conn.size();

		// Do collapsing last, so we can still simplify "trivial loops" like M1-X-L-X-M1
		simplifications += collapseTwoConn(&simplified_net, X_CONN);
	} while(simplifications);
	subtractMols(&mof_asr, &bound_solvent);
	bound_solvent.EndModify();
	mof_asr.EndModify();
	linkers.EndModify();

	if (mil_type_mof) {  // Split 4-coordinated linkers into 3+3 by convention
		if (!fourToTwoThree(&simplified_net, X_CONN)) {
			obErrorLog.ThrowError(__FUNCTION__, "Unexpectedly did not convert 4-coordinated linkers in MIL-like MOF", obWarning);
		}
	}


	// Print out the SMILES for nodes and linkers, and the detected catenation
	analysis << writeFragments(uniqueSMILES(nodes.Separate(), obconv));
	analysis << writeFragments(uniqueSMILES(linkers.Separate(), obconv));
	analysis << "Found " << net_components.size() << " simplified net(s)";

	// Write out the decomposed and simplified MOF, including bond orders
	resetBonds(&nodes);
	resetBonds(&linkers);
	writeCIF(&nodes, "Test/nodes.cif");
	writeCIF(&linkers, "Test/linkers.cif");

	// Write out detected solvents
	writeCIF(&free_solvent, "Test/free_solvent.cif");
	writeCIF(&bound_solvent, "Test/bound_solvent.cif");
	writeCIF(&mof_fsr, "Test/mof_fsr.cif");
	writeCIF(&mof_asr, "Test/mof_asr.cif");

	// Topologically relevant information about the simplified net
	writeCIF(&simplified_net, "Test/removed_two_conn_for_topology.cif");
	writeSystre(&simplified_net, "Test/topology.cgd", X_CONN);

	// Format the fragment keys (psuedo atoms to SMILES) after invalidating unused fragments
	std::map<int,int> active_pseudo_atoms = getNumericFormula(&simplified_net);
	std::vector<int> active_pseudo_elements;
	for (std::map<int,int>::iterator it=active_pseudo_atoms.begin(); it!=active_pseudo_atoms.end(); ++it) {
		active_pseudo_elements.push_back(it->first);
	}
	std::map<std::string,int> removed_keys;
	std::map<std::string,int> gen = node_conv.get_map();
	for (std::map<std::string,int>::iterator it=gen.begin(); it!=gen.end(); ++it) {
		if (!inVector<int>(it->second, active_pseudo_elements)) {
			removed_keys[it->first] = it->second;
		}
	}
	gen = linker_conv.get_map();
	for (std::map<std::string,int>::iterator it=gen.begin(); it!=gen.end(); ++it) {
		if (!inVector<int>(it->second, active_pseudo_elements)) {
			removed_keys[it->first] = it->second;
		}
	}
	for (std::map<std::string,int>::iterator it=removed_keys.begin(); it!=removed_keys.end(); ++it) {
		// Remove unused SMILES codes.  remove_key does not delete anything if the key does not exist
		node_conv.remove_key(it->first);
		linker_conv.remove_key(it->first);
	}
	writeFragmentKeys(node_conv.get_map(), linker_conv.get_map(), removed_keys, X_CONN, "Test/keys_for_condensed_linkers.txt");

	return(analysis.str());
}

extern "C" {
void analyzeMOFc(const char *cifdata, char *analysis, int buflen) {
	// Wrap analyzeMOF with C compatibility for Emscripten usage
	// Make the return value an input param, since c_str is (likely?) a pointer to
	// an internal data structure in std::string, which will go out of scope.
	// A good buflen for the output might be 2^16, or 65536
	const char *TEMP_FILE = "from_emscripten.cif";
	std::ofstream cifp(TEMP_FILE, std::ios::out | std::ios::trunc);
	cifp << std::string(cifdata);
	cifp.close();

	strncpy(analysis, analyzeMOF(TEMP_FILE).c_str(), buflen);
}

int SmilesToSVG(const char* smiles, int options, void* mbuf, unsigned int buflen) {
	// From Noel O'Boyle: https://baoilleach.blogspot.com/2015/02/cheminformaticsjs-open-babel.html
	// See also http://baoilleach.webfactional.com/site_media/blog/emscripten/openbabel/webdepict.html
	OpenBabel::OBMol mol;
	OpenBabel::OBConversion conv;
	conv.SetInFormat("smi");

	bool ok = conv.ReadString(&mol, smiles);
	if (!ok) return -1;

	if (options==0) {
		conv.SetOutFormat("ascii");
		conv.AddOption("a", OpenBabel::OBConversion::OUTOPTIONS, "2.2");

	} else {
		conv.SetOutFormat("svg");
		conv.AddOption("C", OpenBabel::OBConversion::OUTOPTIONS);
		conv.AddOption("P", OpenBabel::OBConversion::OUTOPTIONS, "500");
	}

	std::string out = conv.WriteString(&mol);

	if (out.size()+1 >= buflen)
			return -1;

	char* dst = (char*)mbuf;
	strcpy(dst, out.c_str());
	return out.size();
}
}  // extern "C"


bool readCIF(OBMol* molp, std::string filepath, bool bond_orders, bool makeP1) {
	// Read the first distinguished molecule from a CIF file
	// (TODO: check behavior of mmcif...)
	OBConversion obconversion;
	obconversion.SetInFormat("mmcif");
	obconversion.AddOption("p", OBConversion::INOPTIONS);
	// Defer bond detection until later, once symmetry options are applied
	obconversion.AddOption("b", OBConversion::INOPTIONS);
	bool success = obconversion.ReadFile(molp, filepath);

	if (success && makeP1) {
		// TODO: Consider adding a disorder removal step ("*" and "?" atom labels like CoRE MOF) before applying symmetry operations
		OBUnitCell* uc = (OBUnitCell*)molp->GetData(OBGenericDataType::UnitCell);
		if (!uc) {
			obErrorLog.ThrowError(__FUNCTION__, "Attempted to convert the CIF to P1 without a proper unit cell.", obError);
			success = false;
		} else {
			obErrorLog.ThrowError(__FUNCTION__, "Applying symmetry operations to convert the MOF to P1 (or keep it as P1).", obDebug);
			uc->FillUnitCell(molp);  // since we're operating on a pointer, internal changes to UC should happen automatically
		}
	}

	detectSingleBonds(molp);  // Run single bond detection after filling in the unit cell to avoid running it twice
	detectPaddlewheels(molp);
	if (bond_orders) {
		molp->PerceiveBondOrders();
	}

	return success;
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

OBMol initMOF(OBMol *orig_in_uc) {
	// Initializes a MOF with the same lattice params as *orig_in_uc
	OBMol dest;
	dest.SetData(getPeriodicLattice(orig_in_uc)->Clone(NULL));
	dest.SetPeriodicMol();
	return dest;
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

	OBUnitCell* uc = getPeriodicLattice(pmol);

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
					vector3 unwrap_n_coords = uc->UnwrapFractionalNear(raw_n_coords, x_coords);
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
								vector3 unwrap_v2_coords = uc->UnwrapFractionalNear(raw_v2_coords, unwrap_n_coords);
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
			std::vector<int> bond_dir = GetPeriodicDirection(&*b);
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

void writeFragmentKeys(std::map<std::string,int> nodes, std::map<std::string,int> linkers, std::map<std::string,int> removed, int connector, std::string filepath) {
	// Save fragment identities for condensed_linkers.cif
	std::string equal_line = "================";
	std::ofstream out_file;
	out_file.open(filepath.c_str());
	out_file << "Fragment identities for condensed_linkers.cif" << std::endl << std::endl;

	out_file << "Nodes" << std::endl << equal_line << std::endl;
	for (std::map<std::string,int>::iterator it=nodes.begin(); it!=nodes.end(); ++it) {
		out_file << OBElements::GetSymbol(it->second) << ": " << it->first;
	}
	out_file << std::endl << std::endl;

	out_file << "Linkers" << std::endl << equal_line << std::endl;
	for (std::map<std::string,int>::iterator it=linkers.begin(); it!=linkers.end(); ++it) {
		out_file << OBElements::GetSymbol(it->second) << ": " << it->first;
	}
	out_file << std::endl << std::endl;

	if (X_CONN) {
		out_file << "Connection atom (\"X\")" << std::endl << equal_line << std::endl;
		out_file << OBElements::GetSymbol(X_CONN) << ": " << "<X connector>" << std::endl;
		out_file << std::endl << std::endl;
	}

	if (removed.size()) {
		out_file << "Unused pseudo atom types (see test_partial.cif)" << std::endl << equal_line << std::endl;
		for (std::map<std::string,int>::iterator it=removed.begin(); it!=removed.end(); ++it) {
			out_file << OBElements::GetSymbol(it->second) << ": " << it->first;
		}
		out_file << std::endl << std::endl;
	}

	out_file.close();
}

std::string writeFragments(const std::vector<std::string> &unique_smiles) {
	// Write a list of fragments
	// Use a const_iterator since we're not modifying the vector: http://stackoverflow.com/questions/4890497/how-do-i-iterate-over-a-constant-vector
	// TODO: consider stripping out extraneous tabs, etc, here or elsewhere in the code.
	std::stringstream fragments;
	for (std::vector<std::string>::const_iterator i2 = unique_smiles.begin(); i2 != unique_smiles.end(); ++i2) {
		fragments << *i2;
	}
	return fragments.str();
}

std::string getSMILES(OBMol fragment, OBConversion obconv) {
	// Prints SMILES based on OBConversion parameters
	OBMol canon = fragment;
	resetBonds(&canon);
	unwrapFragmentMol(&canon);
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

void resetBonds(OBMol *mol) {
	// Resets bond orders and bond detection for molecular fragments
	// Starting with a "clean" OBMol is the easiest way to handle this

	std::queue<MinimalAtom> orig_atoms;
	FOR_ATOMS_OF_MOL(a, *mol) {
		MinimalAtom sa;
		sa.loc = a->GetVector();
		sa.element = a->GetAtomicNum();
		sa.is_paddlewheel = a->HasData("Paddlewheel");
		orig_atoms.push(sa);
	}

	OBUnitCell uc_copy = *getPeriodicLattice(mol);
	mol->Clear();
	// Need to allocate memory so it's a new copy that won't go out of scope
	OBUnitCell* uc_data = new OBUnitCell;
	*uc_data = uc_copy;
	mol->SetData(uc_data);
	mol->SetPeriodicMol();

	mol->BeginModify();
	while (!orig_atoms.empty()) {
		MinimalAtom curr_atom = orig_atoms.front();
		orig_atoms.pop();
		OBAtom* copied_atom = formAtom(mol, curr_atom.loc, curr_atom.element);
		if (curr_atom.is_paddlewheel) {
			OBPairData *dp = new OBPairData;
			dp->SetAttribute("Paddlewheel");
			copied_atom->SetData(dp);
		}
		// Consider saving and resetting formal charge as well, e.g. a->SetFormalCharge(0)
	}

	detectSingleBonds(mol);
	// Bond metal atoms in paddlewheels together
	FOR_ATOMS_OF_MOL(a1, *mol) {
		if (a1->HasData("Paddlewheel")) {
			OBAtom* closest_pw = NULL;
			double closest_dist = 100.0;
			FOR_ATOMS_OF_MOL(a2, *mol) {
				if ( a2->HasData("Paddlewheel") && &*a1 != &*a2) {
					double a2_dist = a1->GetDistance(&*a2);
					if (a2_dist < closest_dist) {
						closest_dist = a2_dist;
						closest_pw = &*a2;
					}
				}
			}
			if (!closest_pw) {
				if (mol->NumAtoms() > 1) {  // Do not raise a warning when taking the SMILES of an isolated metal atom
					obErrorLog.ThrowError(__FUNCTION__, "Unable to reconnect a paddlewheel metal to its partner", obWarning);
				}
			} else if (!mol->GetBond(&*a1, closest_pw)) {
				formBond(mol, &*a1, closest_pw);
			}
		}
	}

	mol->PerceiveBondOrders();
	FOR_BONDS_OF_MOL(b, *mol) {
		if ( b->GetBeginAtom()->HasData("Paddlewheel")
			&& b->GetEndAtom()->HasData("Paddlewheel")
			&& b->GetBondOrder() == 3 ) {
			b->SetBondOrder(1);  // Consider normalizing all M#M bonds similarly with isMetal condition
		}
	}

	mol->EndModify();
	normalizeCharges(mol);
}

void detectSingleBonds(OBMol *mol, double skin, bool only_override_oxygen) {
	// Enhances OBMol::ConnectTheDots by also allowing certain cases to exceed maximum valence.
	// By default, the skin (beyond sum of covalent radii) is set to 0.45 AA to match Open Babel.
	// If only_override_oxygen, the only nodular oxygen species get extra valence.  Otherwise, everything.
	// TODO: consider running M-M bonds as another special case besides oxygen.
	// TODO: might also analyze nodes only using using single bonds to avoid related issues with bond order.

	mol->ConnectTheDots();

	const double MIN_DISTANCE = 0.40;
	obErrorLog.ThrowError(__FUNCTION__,
                          "Ran custom detectSingleBonds", obAuditMsg);
	std::vector<OBAtom*> atoms;
	std::vector<double> rads;
	FOR_ATOMS_OF_MOL(a, *mol) {
		atoms.push_back(&*a);
		rads.push_back(OBElements::GetCovalentRad(a->GetAtomicNum()));
	}
	int num_atoms = atoms.size();

	// Note: the N^2 algorithm may be less efficient than Open Babel's version:
	// https://www.slideshare.net/NextMoveSoftware/rdkit-gems
	for (int i = 0; i < num_atoms; ++i) {
		OBAtom* a1 = atoms[i];
		if (only_override_oxygen && a1->GetAtomicNum() != 8) {
			continue;
		}
		std::vector<OBAtom*> nbors_to_bond;
		bool bonded_to_metal = false;

		// In a general neighbor detection algorithm, we would loop j from i+1 to num_atoms.
		// But here, we need to find all neighbors for selected atoms, which is not two-way.
		for (int j = 0; j < num_atoms; ++j) {
			OBAtom* a2 = atoms[j];
			double r = a1->GetDistance(a2);
			double cutoff = rads[i] + rads[j] + skin;

			if (r < cutoff && r > MIN_DISTANCE) {
				if (isMetal(a2)) {
					bonded_to_metal = true;
				}
				if (!(a1->IsConnected(a2))) {
					nbors_to_bond.push_back(a2);
				}
			}
		}

		if (!only_override_oxygen || bonded_to_metal) {
			for (std::vector<OBAtom*>::iterator it = nbors_to_bond.begin(); it != nbors_to_bond.end(); ++it) {
				mol->AddBond(a1->GetIdx(), (*it)->GetIdx(), 1);
			}
		}
	}
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
			obErrorLog.ThrowError(__FUNCTION__, "Attempted to subtract an atom not present in the original molecule", obWarning);
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

bool isSubMol(OBMol *sub, OBMol *super) {
	// Are the atoms in *sub a subset of the atoms in *super?
	FOR_ATOMS_OF_MOL(a, *sub) {
		if (!atomInOtherMol(&*a, super)) {  // NULL if atom not found
			return false;
		}
	}
	return true;
}

std::map<int,int> getNumericFormula(OBMol *mol) {
	// What is the molecular formula, based on atomic numbers?
	// Returns <atomic number of an element, atom count in *mol>
	std::map<int,int> formula;
	FOR_ATOMS_OF_MOL(a, *mol) {
		int element = a->GetAtomicNum();
		if (formula.find(element) == formula.end()) {
			formula[element] = 0;
		}
		formula[element] += 1;
	}
	return formula;
}

ConnExtToInt getLinksToExt(OBMol *mol, OBMol *fragment) {
	// What are the external connections to atoms of a fragment?
	std::vector<OBAtom*> orig_sbu;
	FOR_ATOMS_OF_MOL(a, *fragment) {
		OBAtom* orig_atom = atomInOtherMol(&*a, mol);
		if (!orig_atom) {
			obErrorLog.ThrowError(__FUNCTION__, "Tried to access fragment not present in molecule", obError);
		} else {
			orig_sbu.push_back(orig_atom);
		}
	}

	ConnExtToInt connections;  // <External atom, internal atom bonded to it>
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
	return connections;
}

std::vector<OBAtom*> uniqueExtAtoms(OBMol *mol, OBMol *fragment) {
	// Which unique external atoms does the fragment bind to?
	// Note: this is not necessarily the total number of external connections for a fragment,
	// because it could bind to multiple images of an external atom (e.g. hMOF-0/MOF-5)
	ConnExtToInt connections = getLinksToExt(mol, fragment);
	std::vector<OBAtom*> ext_atoms;
	for (ConnExtToInt::iterator it = connections.begin(); it != connections.end(); ++it) {
		OBAtom* ext = (*it)[0];
		if (!inVector<OBAtom*>(ext, ext_atoms)) {
			ext_atoms.push_back(ext);
		}
	}
	return ext_atoms;
}

MapOfAtomVecs neighborsOverConn(OBAtom *loc, int skip_element) {
	// Gets the nearest neighbors to *loc, skipping over neighbors with an
	// atomic number of skip_element, corresponding to connection site pseudo-atoms
	std::set<OBAtom*> nbors;
	std::set<OBAtom*> skipped;
	std::queue<OBAtom*> to_visit;
	std::set<OBAtom*> visited;

	to_visit.push(loc);
	visited.insert(loc);

	while (!to_visit.empty()) {
		OBAtom* current = to_visit.front();
		to_visit.pop();
		if (current->GetAtomicNum() == skip_element || current == loc) {
			FOR_NBORS_OF_ATOM(n, *current) {
				if (visited.find(&*n) == visited.end()) {
					to_visit.push(&*n);
					visited.insert(&*n);
				}
			}
			if (current->GetAtomicNum() == skip_element) {
				skipped.insert(current);
			}
		} else {
			nbors.insert(current);
		}
	}

	MapOfAtomVecs results;
	results["nbors"] = nbors;
	results["skipped"] = skipped;
	return results;
}

OBAtom* collapseSBU(OBMol *mol, OBMol *fragment, int element, int conn_element) {
	// Simplifies *mol by combining all atoms from *fragment into a single pseudo-atom
	// with atomic number element, maintaining connections to existing atoms.
	// If conn_element != 0, then add that element as a spacer at the connection
	// site (like "X" or "Q" atoms in top-down crystal generators).
	// Returns the pointer to the generated pseudo atom.

	ConnExtToInt connections = getLinksToExt(mol, fragment);

	vector3 centroid = getCentroid(fragment, false);

	mol->BeginModify();

	OBAtom* pseudo_atom = formAtom(mol, centroid, element);
	if (conn_element == 0) {  // Not using an "X" psuedo-atom
		std::vector<OBAtom*> new_external_conn;
		for (ConnExtToInt::iterator it = connections.begin(); it != connections.end(); ++it) {
			OBAtom* external_atom = (*it)[0];
			if (!inVector<OBAtom*>(external_atom, new_external_conn)) {  // Only form external bonds once
				formBond(mol, pseudo_atom, external_atom, 1);
				new_external_conn.push_back(external_atom);
			}
		}
	} else {
		for (ConnExtToInt::iterator it = connections.begin(); it != connections.end(); ++it) {
			// Put the connection point 1/3 of the way between the centroid and the connection midpoint to the exterior
			// (e.g. 1/3 of the way between a BDC centroid and the O-M bond in the -COO group).
			// In a simplified M-X-X-M' system, this will have the convenient property of being mostly equidistant.
			// Note: this follows the convention of many top-down MOF generators placing the connection point halfway on the node-linker bond.
			// In this circumstance, the convention also has the benefit that a linker with many connections to the same metal (-COO)
			// or connections to multiple metals (MOF-74 series) have unique positions for the X_CONN pseudo atoms.
			OBUnitCell* lattice = getPeriodicLattice(mol);
			OBAtom* external_atom = (*it)[0];
			OBAtom* internal_atom = (*it)[1];

			vector3 internal_atom_loc = lattice->UnwrapCartesianNear(internal_atom->GetVector(), centroid);
			vector3 external_atom_loc = lattice->UnwrapCartesianNear(external_atom->GetVector(), internal_atom_loc);
			vector3 conn_loc = lattice->WrapCartesianCoordinate((4.0*centroid + internal_atom_loc + external_atom_loc) / 6.0);

			OBAtom* conn_atom = formAtom(mol, conn_loc, conn_element);
			formBond(mol, conn_atom, pseudo_atom, 1);  // Connect to internal
			formBond(mol, conn_atom, external_atom, 1);
		}
	}

	// Delete SBU from the original molecule
	subtractMols(mol, fragment);

	mol->EndModify();

	return pseudo_atom;
}

int collapseTwoConn(OBMol *net, int ignore_element) {
	// Collapses two-connected nodes into edges to simplify the topology
	// Returns the number of nodes deleted from the network

	net->BeginModify();
	std::vector<OBAtom*> to_delete;
	FOR_ATOMS_OF_MOL(a, *net) {
		if (a->GetValence() == 2 && a->GetAtomicNum() != ignore_element) {
			to_delete.push_back(&*a);
		}
	}

	for (std::vector<OBAtom*>::iterator it = to_delete.begin(); it != to_delete.end(); ++it) {
		std::vector<OBAtom*> nbors;
		FOR_NBORS_OF_ATOM(n, **it) {
			nbors.push_back(&*n);
		}

		if (X_CONN) {  // Transform two-connected node/linker into a connection site
			OBAtom* connector = formAtom(net, (*it)->GetVector(), X_CONN);
			formBond(net, connector, nbors[0]);
			formBond(net, connector, nbors[1]);
		} else {  // Just collapse the node into an edge
			formBond(net, nbors[0], nbors[1]);
		}

		net->DeleteAtom(*it);
	}
	net->EndModify();

	return to_delete.size();
}

int collapseXX(OBMol *net, int element_x) {
	// Simplify X-X bonds in the simplified net with their midpoint
	// Returns the number of X-X bonds simplified
	int simplifications = 0;
	OBUnitCell* lattice = getPeriodicLattice(net);

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
						obErrorLog.ThrowError(__FUNCTION__, nborMsg.str(), obWarning);
					}
					if (x_nbors.size() == 2 && x_nbors[0] == x_nbors[1]) {
						// x1 and x2 are bonded to the same psuedo-atom (M-x1-x2-M').
						// They're not redundant, because M' is often the periodic image of M.
						// Simplifying x1-x2 will cause an inconsistency with two M-X bonds, so skip these x's
						continue;
					}

					// Replace X-X with a new atom at the midpointMake a new atom at the X-X midpoint
					OBAtom* mid_atom = formAtom(net, getMidpoint(x1, x2), element_x);
					for (std::vector<OBAtom*>::iterator it = x_nbors.begin(); it != x_nbors.end(); ++it) {
						formBond(net, mid_atom, *it, 1);
					}
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
	// Remove redundant L-X bonds connecting the same linkers and nodes in the same direction.
	// Simplify the pair of L-X bonds with a new X connector at the midpoint of the two X's.
	// Returns the number of modifications to L-X bonds.

	std::vector<OBAtom*> to_delete;

	FOR_ATOMS_OF_MOL(L, *net) {  // Iterating over linkers
		if (inVector<int>(L->GetAtomicNum(), linker_elements)) {  // if linker atom
			std::vector<OBAtom*> connectors;
			std::map<OBAtom*,OBAtom*> metals;  // <Atom X, M of L-X-M>
			bool found_metals = true;  // Checks if L-X-M format is actually valid
			FOR_NBORS_OF_ATOM(n, *L) {
				if (n->GetAtomicNum() == element_x) {  // connection site
					connectors.push_back(&*n);
					FOR_NBORS_OF_ATOM(M, *n) {  // So we can compare the node connections later
						if (!inVector<int>(M->GetAtomicNum(), linker_elements) && M->GetAtomicNum() != element_x) {
							metals[&*n] = &*M;
						}
					}
					if (metals.find(&*n) == metals.end()) {
						found_metals = false;
					}
				}
			}
			if (!found_metals) {
				obErrorLog.ThrowError(__FUNCTION__, "Cannot simplify LX: no metal at other end of L-X-M", obDebug);
				continue;
			}

			// Iterate through the bonded L-X's to find redundant pairs: two X's
			// that connect to the same node in the same unit cell are deleted here.
			for (std::vector<OBAtom*>::iterator x1 = connectors.begin(); x1 != connectors.end(); ++x1) {
				for (std::vector<OBAtom*>::iterator x2 = connectors.begin(); x2 != connectors.end(); ++x2) {
					// By merit of the FOR loop over L, we've already established that L is the same.  Also check M with the metals map.
					if ( *x1 != *x2
						&& !inVector<OBAtom*>(*x1, to_delete)
						&& !inVector<OBAtom*>(*x2, to_delete)
						&& metals[*x1] == metals[*x2] )
					{
						OBMol x_test = initMOF(net);
						x_test.BeginModify();
						std::vector<OBAtom*> test_atoms;  // Copy of M1-X1-L-X2-M2, where M2 might equal M1
						test_atoms.push_back(formAtom(&x_test, metals[*x1]->GetVector(), metals[*x1]->GetAtomicNum()));
						test_atoms.push_back(formAtom(&x_test, (*x1)->GetVector(), (*x1)->GetAtomicNum()));
						test_atoms.push_back(formAtom(&x_test, (&*L)->GetVector(), (&*L)->GetAtomicNum()));
						test_atoms.push_back(formAtom(&x_test, (*x2)->GetVector(), (*x2)->GetAtomicNum()));
						// M2 atom is already defined by metals[*x1].  See parent conditional statement.
						formBond(&x_test, test_atoms[0], test_atoms[1], 1);
						formBond(&x_test, test_atoms[1], test_atoms[2], 1);
						formBond(&x_test, test_atoms[2], test_atoms[3], 1);
						formBond(&x_test, test_atoms[3], test_atoms[0], 1);
						x_test.EndModify();

						// If x_test is periodic, then M2 is in a different UC than M1, thus X1-X2 is a bridge.
						if (!isPeriodicChain(&x_test)) {
							to_delete.push_back(*x1);
							to_delete.push_back(*x2);
							OBAtom* x_mid = formAtom(net, getMidpoint(*x1, *x2, false), element_x);
							formBond(net, x_mid, metals[*x1], 1);
							formBond(net, x_mid, &*L, 1);
						}
					}
				}
			}
		}
	}

	// Delete the redundant (original) X's (which will also remove their X-L and X-M bonds)
	net->BeginModify();
	for (std::vector<OBAtom*>::iterator it = to_delete.begin(); it != to_delete.end(); ++it) {
		net->DeleteAtom(*it);
	}
	net->EndModify();

	return to_delete.size();
}

std::vector<OBMol> removeOneConn(OBMol *net, std::map<OBAtom*, OBMol*> pseudo_to_real, std::vector<int> allowed_elements, int connector) {
	// Remove ligands one-connected to a node.
	// Currently classify these as bound solvents, but they could represent
	// post-synthetic modification or organic parts of the node.
	// Assumes that MOF nodes are always metal oxides.
	// Inputs are the simplified net, a conversion from net pseudo atoms to molecules of real atoms,
	// the elements allowed for removal, and the connector atom type.
	// Returns the vector of OBMol's deleted from the simplified net, and cleans up pseudo_to_real map.

	std::vector<OBMol> fragments_removed;
	std::vector<OBAtom*> pseudo_deletions;
	net->BeginModify();

	FOR_ATOMS_OF_MOL(a, *net) {
		if (a->GetValence() == 1) {
			int point_type = a->GetAtomicNum();
			if (inVector<int>(point_type, allowed_elements) || allowed_elements.size() == 0) {  // ignore check if not specified
				if (X_CONN) { // Remove adjacent connection points, simplifying the problem to the base case
					std::set<OBAtom*> a_nbors = neighborsOverConn(&*a, X_CONN)["nbors"];
					std::set<OBAtom*> a_conns = neighborsOverConn(&*a, X_CONN)["skipped"];
					for (std::set<OBAtom*>::iterator it=a_conns.begin(); it!=a_conns.end(); ++it) {
						pseudo_deletions.push_back(*it);  // Delete the relevant connector atoms at the end
					}
					formBond(net, *(a_nbors.begin()), &*a);
				}

				// Save the bound ligand.  At the end, we'll delete the one-connected "nodes"
				fragments_removed.push_back(*(pseudo_to_real[&*a]));
				pseudo_deletions.push_back(&*a);
				pseudo_to_real.erase(&*a);  // Map entry will be invalid after pseudo atom is deleted
			} else {
				obErrorLog.ThrowError(__FUNCTION__, "Unexpected one-connected component after net simplifications.", obWarning);
			}
		}
	}

	for (std::vector<OBAtom*>::iterator it=pseudo_deletions.begin(); it!=pseudo_deletions.end(); ++it) {
		net->DeleteAtom(*it);
	}
	net->EndModify();

	return fragments_removed;
}

int fourToTwoThree(OBMol *net, int X_CONN) {
	// Transforms four-connected atoms to two three-connected pseudo atoms.
	// This satisfies the convention used for MIL-47-like topologies.
	// Returns the number of simplified linkers

	int changes = 0;
	std::queue<OBAtom*> to_change;

	FOR_ATOMS_OF_MOL(a, *net) {
		if (a->GetValence() == 4) {
			to_change.push(&*a);
		}
	}

	net->BeginModify();
	while (!to_change.empty()) {
		OBAtom* current = to_change.front();
		to_change.pop();
		std::vector<OBAtom*> nbors;
		FOR_NBORS_OF_ATOM(n, *current) {
			nbors.push_back(&*n);
		}

		// The four-connected node will have two sides to break apart the one node into two
		std::vector<OBAtom*> side1;
		side1.push_back(nbors[0]);
		side1.push_back(minAngleNbor(current, side1[0]));
		std::vector<OBAtom*> side2;
		for(std::vector<OBAtom*>::iterator it=nbors.begin(); it!=nbors.end(); ++it) {
			if (!inVector<OBAtom*>(*it, side1)) {
				side2.push_back(*it);
			}
		}

		if ((side1[0])->GetAngle(current, side1[1]) > 85) {
			obErrorLog.ThrowError(__FUNCTION__, "Trying to simplify a square-like four-connected pseudo atom", obWarning);
		}
		bool mismatch = false;  // Verify that the pairings are self-consistent for all four atoms
		for (std::vector<OBAtom*>::iterator it=nbors.begin(); it!=nbors.end(); ++it) {
			if (inVector<OBAtom*>(*it, side1)) {
				if (!inVector<OBAtom*>(minAngleNbor(current, *it), side1)) {
					mismatch = true;
				}
			} else {
				if (inVector<OBAtom*>(minAngleNbor(current, *it), side1)) {
					mismatch = true;
				}
			}
		}
		if (mismatch) {
			obErrorLog.ThrowError(__FUNCTION__, "Mismatched neighbors when assigning the split", obWarning);
		}

		int pseudo_element = current->GetAtomicNum();
		vector3 current_loc = current->GetVector();
		vector3 pseudo_loc;
		OBAtom* pseudo_1 = formAtom(net, getMidpoint(side1[0], side1[1]), pseudo_element);
		OBAtom* pseudo_2 = formAtom(net, getMidpoint(side2[0], side2[1]), pseudo_element);
		formBond(net, pseudo_1, side1[0], 1);
		formBond(net, pseudo_1, side1[1], 1);
		formBond(net, pseudo_2, side2[0], 1);
		formBond(net, pseudo_2, side2[1], 1);

		if (X_CONN) {
			OBAtom* conn_atom = formAtom(net, current_loc, X_CONN);
			formBond(net, pseudo_1, conn_atom, 1);
			formBond(net, pseudo_2, conn_atom, 1);
		} else {  // Note: usage without X_CONN is untested.
			formBond(net, pseudo_1, pseudo_2, 1);
		}
		net->DeleteAtom(current);
		changes += 1;
	}
	net->EndModify();

	return changes;
}

OBAtom* minAngleNbor(OBAtom* base, OBAtom* first_conn) {
	// Which neighbor to base has the smallest first-base-nbor bond angle?
	// (excluding first_conn, of course)
	OBAtom* min_nbor = NULL;
	double min_angle = 360.0;
	FOR_NBORS_OF_ATOM(n, *base) {
		if (&*n != first_conn) {
			double test_angle = first_conn->GetAngle(base, &*n);
			if (test_angle < min_angle) {
				min_angle = test_angle;
				min_nbor = &*n;
			}
		}
	}
	return min_nbor;
}

UCMap unwrapFragmentUC(OBMol *fragment, bool allow_rod, bool warn_rod) {
	// Starting with a random atom, traverse the fragment atom-by-atom to determine
	// which unit cell each atom belongs to in a self-consistent manner.
	// Includes an optional parameter to allow 1D periodic fragments (e.g. MIL-47).
	// By default, these are forbidden (since the ordering is undefined) and returns an empty map.

	std::queue<OBAtom*> to_visit;
	std::map<OBAtom*, std::vector<int> > unit_cells;
	// Start at whichever atom is (randomly?) saved first
	// Note: atom arrays begin with 1 in OpenBabel, while bond arrays begin with 0.
	OBAtom* start_atom = fragment->GetAtom(1);
	to_visit.push(start_atom);
	unit_cells[start_atom] = makeVector(0, 0, 0);  // original unit cell

	while (!to_visit.empty()) {
		OBAtom* current = to_visit.front();
		to_visit.pop();
		FOR_NBORS_OF_ATOM(nbr, current) {
			OBBond* nbr_bond = fragment->GetBond(current, &*nbr);
			std::vector<int> uc = GetPeriodicDirection(nbr_bond);  // TODO: Is there a good way to remove this function?  Otherwise, add it as a MOF helper method
			if (nbr_bond->GetBeginAtom() == &*nbr) {  // opposite bond direction as expected
				uc = makeVector(-1*uc[0], -1*uc[1], -1*uc[2]);
			}

			std::vector<int> current_uc = unit_cells[current];
			uc = makeVector(current_uc[0] + uc[0], current_uc[1] + uc[1], current_uc[2] + uc[2]);

			if (unit_cells.find(&*nbr) == unit_cells.end()) {  // Unvisited atom
				// Make sure to visit the neighbor (and its neighbors, etc.)
				// Each atom will only be traversed once, since we've already added it to unit_cells
				to_visit.push(&*nbr);
				unit_cells[&*nbr] = uc;
			} else {  // Visited atom: check for loops across periodic boundaries
				if (unit_cells[&*nbr] != uc) {
					if (warn_rod) {
						obErrorLog.ThrowError(__FUNCTION__, "Found periodic loop when unwrapping fragment.  Unit cells are may not be self-consistent.", obWarning);
					}
					if (!allow_rod) {
						unit_cells.clear();
						return unit_cells;
					}
				}
			}
		}
	}

	if (fragment->NumAtoms() != unit_cells.size()) {  // Note: will not run if periodic loops are found and !allow_rod
		obErrorLog.ThrowError(__FUNCTION__, "More than one fragment found.  Behavior is undefined.", obError);
		unit_cells.clear();
	}
	return unit_cells;
}

bool unwrapFragmentMol(OBMol* fragment) {
	// Starting with a random atom in a fragment, unwrap the atomic coordinates
	// to all belong to the same unit cell.
	// Modifies fragment unless it contains a rod; otherwise return false.
	// TODO: consider refactoring getCentroid and other codes related to UCMap?

	UCMap rel_uc = unwrapFragmentUC(fragment, false, false);
	if (rel_uc.size() == 0) {
		return false;
	}

	for (UCMap::iterator it=rel_uc.begin(); it!=rel_uc.end(); ++it) {
		OBAtom* curr_atom = it->first;
		std::vector<int> uc_shift = it->second;

		vector3 uc_shift_frac(uc_shift[0], uc_shift[1], uc_shift[2]);  // Convert ints to doubles
		vector3 coord_shift = getPeriodicLattice(fragment)->FractionalToCartesian(uc_shift_frac);
		curr_atom->SetVector(curr_atom->GetVector() + coord_shift);
	}
	return true;
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

	// The more complicated periodic case requires "unwrapping" the molecular fragment
	UCMap unit_cells = unwrapFragmentUC(fragment, true, true);
	OBUnitCell* lattice = getPeriodicLattice(fragment);
	for (UCMap::iterator it=unit_cells.begin(); it!=unit_cells.end(); ++it) {
		double weight = 1.0;
		if (weighted) {
			weight = it->first->GetAtomicMass();
		}
		std::vector<int> uc_shift = it->second;  // <first: second> = <key: value> of a map/dict.
		vector3 uc_shift_frac(uc_shift[0], uc_shift[1], uc_shift[2]);  // Convert ints to doubles
		vector3 coord_shift = lattice->FractionalToCartesian(uc_shift_frac);
		center += weight * (it->first->GetVector() + coord_shift);
		total_weight += weight;
	}
	center /= total_weight;
	center = lattice->WrapCartesianCoordinate(center);  // Keep result in the [0,1] unit cell
	return center;
}

vector3 getMidpoint(OBAtom* a1, OBAtom* a2, bool weighted) {
	// Returns the midpoint between two atoms, accounting for PBC if present

	vector3 a1_raw = a1->GetVector();
	vector3 a2_raw = a2->GetVector();
	double wt_1 = 1.0;
	double wt_2 = 1.0;
	if (weighted) {
		wt_1 = a1->GetAtomicMass();
		wt_2 = a2->GetAtomicMass();
	}

	OBUnitCell* lattice = getPeriodicLattice(a1->GetParent());
	if (lattice) {
		vector3 a2_unwrap = lattice->UnwrapCartesianNear(a2_raw, a1_raw);
		return lattice->WrapCartesianCoordinate((wt_1 * a1_raw + wt_2 * a2_unwrap) / (wt_1 + wt_2));
	} else {
		return (wt_1 * a1_raw + wt_2 * a2_raw) / (wt_1 + wt_2);
	}
}

bool isPeriodicChain(OBMol *mol) {
	// Is this molecule/fragment 1D periodic?
	// For example, 1D periodic rods (MIL-47 and related topologies).
	// Implements an algorithm similar to dimensionality of channel systems in Zeo++.
	// See description in Sections 2.2.2-2.2.3 of 10.1016/j.micromeso.2011.08.020

	if (!mol || mol->NumAtoms() == 0) {
		return false;
	}
	if (mol->Separate().size() != 1) {
		obErrorLog.ThrowError(__FUNCTION__, "Not processing multi-fragment OBMol", obWarning);
		return false;
	}
	// Unwrapping the fragment will return an empty map if there's UC inconsistencies,
	// which are indicative of a periodic chain (or pore).
	if (unwrapFragmentUC(mol, false, false).size() == 0) {
		return true;
	}
	return false;
}

int sepPeriodicChains(OBMol *nodes) {
	// Separate periodic chains (rods like MIL-47) by disconnecting the M1-O-M2 bonds in *nodes.
	// Returns the number of rods detected and simplified.

	int simplifications = 0;
	std::vector<OBBond*> delbonds;

	std::vector<OBMol> sep_nodes = nodes->Separate();
	for (std::vector<OBMol>::iterator it = sep_nodes.begin(); it != sep_nodes.end(); ++it) {
		if (isPeriodicChain(&*it)) {
			FOR_ATOMS_OF_MOL(a, *it) {
				if (a->GetAtomicNum() == 8) {
					bool all_metals = true;
					FOR_NBORS_OF_ATOM(n, *a) {
						if (!isMetal(&*n)) {
							all_metals = false;
						}
					}
					if (all_metals) {
						OBAtom* orig_atom = atomInOtherMol(&*a, nodes);
						FOR_BONDS_OF_ATOM(b, *orig_atom) {
							if (!inVector(&*b, delbonds)) {
								delbonds.push_back(&*b);
							}
						}
					}
				}
			}
			simplifications += 1;
		}
	}

	nodes->BeginModify();
	for (std::vector<OBBond*>::iterator it=delbonds.begin(); it!=delbonds.end(); ++it) {
		nodes->DeleteBond(*it);
	}
	nodes->EndModify();
	if (simplifications) {  // Check that the periodic nodes are actually separated
		std::vector<OBMol> fixed_sep_nodes = nodes->Separate();
		for (std::vector<OBMol>::iterator it=fixed_sep_nodes.begin(); it!=fixed_sep_nodes.end(); ++it) {
			if (isPeriodicChain(&*it)) {
				std::string err_msg = "Node with formula " + it->GetFormula() + " still periodic after separations";
				obErrorLog.ThrowError(__FUNCTION__, err_msg, obWarning);
			}
		}
	}
	return simplifications;
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

bool normalizeCharges(OBMol *mol) {
	// Correct formal charges on carboxylic acids, imidazolate, etc.
	// Returns if any changes were made to the molecule.

	bool changed = false;
	std::queue<std::pair<std::string, std::string> > reactions;
	// Carboxylate: "OD1-0" means oxygen with one explicit connection, zero charge
	reactions.push(std::make_pair("O=C[OD1-0:1]", "O=C[O-:1]"));
	// FIXME: imidazolate is still broken due to presence of radicals.
	// Newer upstream work on aromaticity/spin handling may be necessary
	// reactions.push(std::make_pair("c1ncc[n:1]1", "c1ncc[n-:1]1"));  // Imidazolate, etc.
	reactions.push(std::make_pair("[nD2:1]1nc[cD3]c1", "[n-:1]1nc[cD3]c1"));  // pyrazole, e.g. sym_3_mc_0

	// TODO: While aromaticity bugs are being sorted out (esp. with nitrogens),
	// we could also consider hard-coding the known nitrogen rings.
	// https://en.wikipedia.org/wiki/Heterocyclic_compound is a good resource.

	while (!reactions.empty()) {
		std::string reactants = reactions.front().first;
		std::string products = reactions.front().second;
		reactions.pop();

		OBChemTsfm tsfm;
		if (!tsfm.Init(reactants, products)) {
			obErrorLog.ThrowError(__FUNCTION__, "Internal error: could not parse reaction transform", obError);
			return false;
		}
		if (tsfm.Apply(*mol)) {
			changed = true;
		}
	}

	return changed;
}

bool detectPaddlewheels(OBMol *mol) {
	// Normalize all paddlewheel bonds to be a single bond between metals, regardless of exact distance.
	// Also sets a "Paddlewheel" attribute on the relevant atoms for reperception of the relevant bond.
	// See also earlier tests and example code from: http://openbabel.org/dev-api/classOpenBabel_1_1OBSmartsPattern.shtml
	// Returns if any paddlewheels were detected in the structure

	bool found_pw = false;
	OBSmartsPattern paddlewheel;
	// Paddlewheel metals might have a M-M bond, and possibly a coordinated solvent or pillar linker
	paddlewheel.Init("[D4,D5,D6:1](OCO1)(OCO2)(OCO3)OCO[D4,D5,D6:2]123");
	paddlewheel.Match(*mol);
	std::vector<std::vector<int> > maplist = paddlewheel.GetUMapList();

	std::vector<std::vector<int> >::iterator i;
	std::vector<int>::iterator j;
	for (i=maplist.begin(); i!=maplist.end(); ++i) {  // loop over matches
		OBMol candidate = initMOF(mol);  // Have to check the match for infinite rods
		candidate.BeginModify();
		std::vector<OBAtom*> pw_metals;
		for (j=i->begin(); j!=i->end(); ++j) {  // loop over paddlewheel atoms
			OBAtom* curr_atom = mol->GetAtom(*j);
			if (isMetal(curr_atom)) {
				pw_metals.push_back(curr_atom);
			}
			formAtom(&candidate, curr_atom->GetVector(), curr_atom->GetAtomicNum());
		}
		candidate.EndModify();
		resetBonds(&candidate);

		// Delete OOC=COO bonds, which would only be present for adjacent paddlewheels
		OBSmartsPattern adjacent_pw;
		adjacent_pw.Init("OC(O)=C(O)O");
		if (adjacent_pw.Match(candidate)) {
			std::vector<std::vector<int> > adjacent_maplist = adjacent_pw.GetUMapList();
			std::vector<OBBond*> adj_delete;
			std::vector<std::vector<int> >::iterator match_adj;

			for (match_adj=adjacent_maplist.begin(); match_adj!=adjacent_maplist.end(); ++match_adj) {
				OBAtom* c1 = candidate.GetAtom((*match_adj)[1]);
				OBAtom* c2 = candidate.GetAtom((*match_adj)[3]);
				if (c1->GetAtomicNum() != 6 || c2->GetAtomicNum() != 6) {
					obErrorLog.ThrowError(__FUNCTION__, "Internal error: unexpected atom type matched in SMARTS pattern", obError);
				}
				adj_delete.push_back(candidate.GetBond(c1, c2));
			}

			candidate.BeginModify();
			for (std::vector<OBBond*>::iterator it=adj_delete.begin(); it!=adj_delete.end(); ++it) {
				obErrorLog.ThrowError(__FUNCTION__, "Deleted bond between adjacent paddlewheel carboxylates in candidate match", obDebug);
				candidate.DeleteBond(*it);
			}
			candidate.EndModify();
		}

		if (pw_metals.size() != 2) {
			obErrorLog.ThrowError(__FUNCTION__, "Inconsistent paddlewheel match without two metal atoms", obError);
			return false;
		}

		if (isPeriodicChain(&candidate)) {
			obErrorLog.ThrowError(__FUNCTION__, "Skipping paddlewheel assignment: match is an infinite rod", obDebug);
		} else {
			obErrorLog.ThrowError(__FUNCTION__, "Found a paddlewheel.  Assigining \"Paddlewheel\" attribute", obDebug);
			found_pw = true;
			OBBond* old_bond = mol->GetBond(pw_metals[0], pw_metals[1]);
			if (old_bond) {
				mol->DeleteBond(old_bond);
			}
			formBond(mol, pw_metals[0], pw_metals[1], 1);

			OBPairData *dp = new OBPairData;
			dp->SetAttribute("Paddlewheel");
			//dp->SetValue("some value");
			//dp->SetOrigin(external);
			pw_metals[0]->SetData(dp);
			dp = new OBPairData;
			dp->SetAttribute("Paddlewheel");
			pw_metals[1]->SetData(dp);
		}
	}
	return found_pw;
}

std::vector<int> GetPeriodicDirection(OBBond *bond) {
	// What is the unit cell of the end atom wrt the first?
	// Returns {0,0,0} if not periodic or if wrapping is not required.

	std::vector<int> direction;
	direction.push_back(0);
	direction.push_back(0);
	direction.push_back(0);

	if (bond->IsPeriodic())  // Otherwise, return all zeros
	{
		OBUnitCell *box = getPeriodicLattice(bond->GetParent());
		vector3 begin, end_orig, end_expected, uc_direction;
		begin = box->CartesianToFractional(bond->GetBeginAtom()->GetVector());
		end_orig = box->CartesianToFractional(bond->GetEndAtom()->GetVector());
		end_expected = box->UnwrapFractionalNear(end_orig, begin);

		// To get the signs right, consider the example {0, 0.7}.  We want -1 as the periodic direction.
		// TODO: Think about edge cases, particularly atoms on the border of the unit cell.
		uc_direction = end_expected - end_orig;

		for (int i = 0; i < 3; ++i) {
			double raw_cell = uc_direction[i];
			direction[i] = static_cast<int>(lrint(raw_cell));
		}
	}
	return direction;
}

OBUnitCell* getPeriodicLattice(OBMol *mol) {
	// Replacement for the old OBMol.GetPeriodicLattice helper function
	return (OBUnitCell*)mol->GetData(OBGenericDataType::UnitCell);
}
