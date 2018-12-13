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
#include <openbabel/elements.h>

#include "config_sbu.h"
#include "invector.h"
#include "obdetails.h"
#include "periodic.h"
#include "framework.h"

using namespace OpenBabel;  // See http://openbabel.org/dev-api/namespaceOpenBabel.shtml


// Connections to an atom or fragment, in the form of a set where each element
// is the vector of OBAtom pointers: <External atom, internal atom bonded to it>
typedef std::set<std::vector<OBAtom*> > ConnExtToInt;

typedef std::map<std::string,std::set<OBAtom*> > MapOfAtomVecs;


// Function prototypes
std::string analyzeMOF(std::string filename);
extern "C" void analyzeMOFc(const char *cifdata, char *analysis, int buflen);
extern "C" int SmilesToSVG(const char* smiles, int options, void* mbuf, unsigned int buflen);

void writeSystre(OBMol* molp, std::string filepath, int element_x = 0, bool write_centers = true);
void writeFragmentKeys(std::map<std::string,int> nodes, std::map<std::string,int> linkers, std::map<std::string,int> removed, int connector, std::string filepath);
std::string writeFragments(std::vector<OBMol> fragments, OBConversion obconv);
std::string getSMILES(OBMol fragment, OBConversion obconv);
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
int sepPeriodicChains(OBMol *nodes);


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
	if (!importCIF(&orig_mol, filename, false)) {
		std::cerr << "Error reading file: %s" << filename << std::endl;
		return "";
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
	OBMol free_solvent = initMOFwithUC(&orig_mol);
	OBMol bound_solvent = initMOFwithUC(&orig_mol);

	// TODO: Write out the original mol like this (near instantly) for debugging (maybe as part of importCIF)
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
	// TODO: this is where I could add other branch point detection, such as phenyl rings


	// Print out the SMILES for nodes and linkers, and the detected catenation
	analysis << writeFragments(nodes.Separate(), obconv);
	analysis << writeFragments(linkers.Separate(), obconv);
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
			int3 bond_dir = GetPeriodicDirection(&*b);
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

std::string writeFragments(std::vector<OBMol> fragments, OBConversion obconv) {
	// Write a list of unique SMILES for a set of fragments
	// TODO: consider stripping out extraneous tabs, etc, here or elsewhere in the code.
	std::stringstream written;
	std::set<std::string> unique_smiles;
	for (std::vector<OBMol>::iterator it = fragments.begin(); it != fragments.end(); ++it) {
		unique_smiles.insert(getSMILES(*it, obconv));  // only adds unique values in a set
	}
	for (std::set<std::string>::iterator i2 = unique_smiles.begin(); i2 != unique_smiles.end(); ++i2) {
		written << *i2;
	}

	return written.str();
}

std::string getSMILES(OBMol fragment, OBConversion obconv) {
	// Prints SMILES based on OBConversion parameters
	OBMol canon = fragment;
	resetBonds(&canon);
	unwrapFragmentMol(&canon);
	return obconv.WriteString(&canon);
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
						OBMol x_test = initMOFwithUC(net);
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

