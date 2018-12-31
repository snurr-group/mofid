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
#include "framework.h"
#include "periodic.h"
#include "pseudo_atom.h"
#include "topology.h"
#include "virtual_mol.h"


using namespace OpenBabel;  // See http://openbabel.org/dev-api/namespaceOpenBabel.shtml


typedef std::map<std::string,std::set<OBAtom*> > MapOfAtomVecs;


// Function prototypes
std::string analyzeMOF(std::string filename);
extern "C" void analyzeMOFc(const char *cifdata, char *analysis, int buflen);
extern "C" int SmilesToSVG(const char* smiles, int options, void* mbuf, unsigned int buflen);

void writeFragmentKeys(std::map<std::string,int> nodes, std::map<std::string,int> linkers, std::map<std::string,int> removed, int connector, std::string filepath);
std::string writeFragments(std::vector<OBMol> fragments, OBConversion obconv);
std::string getSMILES(OBMol fragment, OBConversion obconv);
int collapseTwoConn(OBMol* net, int ignore_element = 0);
std::vector<OBMol> removeOneConn(OBMol *net, std::map<OBAtom*, OBMol*> pseudo_to_real, std::vector<int> allowed_elements, int connector);
// Replaced methods below. TODO remove:
// Connections to an atom or fragment, in the form of a set where each element
// is the vector of OBAtom pointers: <External atom, internal atom bonded to it>
typedef std::set<std::vector<OBAtom*> > ConnExtToInt;
ConnExtToInt getLinksToExt(OBMol *mol, OBMol *fragment);
std::vector<OBAtom*> uniqueExtAtoms(OBMol *mol, OBMol *fragment);
MapOfAtomVecs neighborsOverConn(OBAtom *loc, int skip_element);  // entirely deprecated, as well as collapseXX, if we don't have extraneous bonds in our network
// these changes to XX will also greatly simplify Systre

/* Define global parameters for MOF decomposition */
// Atom type for connection sites.  Assigned to Te (52) for now.  Set to zero to disable.
const int X_CONN = 52;

/*
// TODO: delete after implementing the code in topology.cpp
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
*/

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
	/* TODO: rewrite most of this code in Deconstructor and derived classes, saving "main"
	for parsing input args.  MOFidDeconstructor, AllNodeDeconstructor, SingleNodeDeconstructor */

	std::stringstream analysis;
	OBMol orig_mol;
	// Massively improving performance by skipping kekulization of the full MOF
	if (!importCIF(&orig_mol, filename, false)) {
		std::cerr << "Error reading file: %s" << filename << std::endl;
		return "";
	}

	Topology simplified(&orig_mol);
	OBMol testing_output_mol = simplified.ToOBMol();
	writeCIF(&testing_output_mol, "Test/simplified_test_orig.cif");
	// TODO: implement FSR and ASR later
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

	OBMol split_mol;  // temporary molecule to setup the initial metal-based fragmentation
	copyMOF(&orig_mol, &split_mol);

	// TODO: does topology.h, etc., use Begin/EndModify routines?

	// Find linkers by deleting bonds to metals
	// TODO: in this block, for SBU decomposition algorithms, do some manipulations to modify/restore bonds before fragmentation.
	// That will probably take the form of an optional preprocessing step before fragment assignment.
	// Similarly, there will probably be a fragmenter that breaks apart the nodes/linkers using a standard algorithm for node/linker SMILES names.
	std::vector<OBMol> fragments;
	deleteBonds(&split_mol, true);
	fragments = split_mol.Separate();

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
	std::stringstream nonmetalMsg;
	for (std::vector<OBMol>::iterator it = fragments.begin(); it != fragments.end(); ++it) {
		std::string mol_smiles = getSMILES(*it, obconv);
		VirtualMol fragment_act_atoms(&orig_mol);
		fragment_act_atoms.ImportCopiedFragment(&*it);
		VirtualMol fragment_pa = simplified.OrigToPseudo(fragment_act_atoms);;

		// If str comparisons are required, include a "\t\n" in the proposed smiles
		bool all_oxygens = true;  // Also allow hydroxyls, etc.
		FOR_ATOMS_OF_MOL(a, *it){
			if (a->GetAtomicNum() != 8 && a->GetAtomicNum() != 1) {
				all_oxygens = false;
			}
		}

		if (fragment_act_atoms.GetExternalBonds().size() == 0) {
			// Assume free solvents are organic (or lone metals), so they'd be isolated without any external connections
			nonmetalMsg << "Deleting free solvent " << mol_smiles;
			free_solvent += *it;  // TODO: just delete in the simplified object and assign the requisite roles
			// FIXME: revisit this loop and implementation in the main loop as well
			AtomSet free_set = fragment_pa.GetAtoms();
			for (AtomSet::iterator free_it=free_set.begin(); free_it!=free_set.end(); ++free_it) {
				simplified.DeleteAtomAndConns(*free_it);
			}
			// TODO: simplify the net as well
			// (and think about improving or removing Topology::RemoveOrigAtoms)
			// Let's get the general test cases finished first, then come back and implement this one.
		} else if (it->NumAtoms() == 1) {
			nonmetalMsg << "Found a solitary atom with atomic number " << it->GetFirstAtom()->GetAtomicNum() << std::endl;
			simplified.SetRoleToAtoms( "node", fragment_pa);
		} else if (all_oxygens) {
			nonmetalMsg << "Found an oxygen species " << mol_smiles;
			simplified.SetRoleToAtoms("node", fragment_pa);
			// do we condense these yet?  probably later
		} else {
			nonmetalMsg << "Deleting linker " << mol_smiles;
			PseudoAtom collapsed = simplified.CollapseFragment(fragment_pa);
			simplified.SetRoleToAtom("linker", collapsed);
		}
	}
	obErrorLog.ThrowError(__FUNCTION__, nonmetalMsg.str(), obDebug);
	OBMol test_partial = simplified.ToOBMol();
	writeCIF(&test_partial, "Test/test_partial.cif");

	// Simplify all the node SBUs into single points.
	bool mil_type_mof = false;
	VirtualMol node_pa = simplified.GetAtomsOfRole("node");
	node_pa = simplified.FragmentWithIntConns(node_pa);
	std::vector<VirtualMol> node_fragments = node_pa.Separate();
	for (std::vector<VirtualMol>::iterator it=node_fragments.begin(); it!=node_fragments.end(); ++it) {
		VirtualMol fragment_mol = *it;

		OBMol fragment_obmol = fragment_mol.ToOBMol();
		if (!isPeriodicChain(&fragment_obmol)) {  // normal, nonperiodic case
			PseudoAtom collapsed = simplified.CollapseFragment(*it);
			simplified.SetRoleToAtom("node", collapsed);
		} else {  // based on sepPeriodicChains
			// TODO: consider refactoring this code to a method within topology.cpp
			// or the upcoming simplification class
			obErrorLog.ThrowError(__FUNCTION__, "Detecting infinite chains", obInfo);
			mil_type_mof = true;

			// Detect single-atom nonmetal bridging atoms
			AtomSet bridging_atoms;
			AtomSet rod_atoms = it->GetAtoms();
			for (AtomSet::iterator rod_it=rod_atoms.begin(); rod_it!=rod_atoms.end(); ++rod_it) {
				AtomSet rod_orig = simplified.PseudoToOrig(VirtualMol(*rod_it)).GetAtoms();
				if (rod_orig.size() == 1) {
					OBAtom* single_atom = *(rod_orig.begin());
					if (!isMetal(single_atom)) {
						bridging_atoms.insert(*rod_it);  // the PA, not single_atom from the original MOF
					}
				}
			}

			// Handle bridging atoms separately from the rest of the rod.
			for (AtomSet::iterator br_it=bridging_atoms.begin(); br_it!=bridging_atoms.end(); ++br_it) {
				fragment_mol.RemoveAtom(*br_it);
				simplified.SetRoleToAtom("node bridge", *br_it);
				simplified.SetRoleToAtom("node", *br_it, false);
			}
			fragment_mol = simplified.FragmentWithoutConns(fragment_mol);

			if (fragment_mol.NumAtoms() == 0) {
				obErrorLog.ThrowError(__FUNCTION__, "Unexpectedly deleted all atoms in a periodic rod during simplificaiton.", obError);
				continue;
			}

			// Simplify the non-bridging metals
			std::vector<VirtualMol> rod_fragments = fragment_mol.Separate();
			for (std::vector<VirtualMol>::iterator frag_it=rod_fragments.begin(); frag_it!=rod_fragments.end(); ++frag_it) {
				if (frag_it->NumAtoms() > 1) {
					obErrorLog.ThrowError(__FUNCTION__, "Combining metal atoms within a periodic rod (likely okay, but untested code--check it).", obError);
					fragment_mol = simplified.FragmentWithoutConns(fragment_mol);
					PseudoAtom collapsed = simplified.CollapseFragment(fragment_mol);
					simplified.SetRoleToAtom("node", collapsed);
				}  // else, if a single-metal fragment, there's nothing to simplify

			}
		}
	}

	OBMol test_nodes = simplified.ToOBMol();
	writeCIF(&test_nodes, "Test/test_with_simplified_nodes.cif");

	// Simplify the topological net
	int simplifications = 0;
	do {
		simplifications = 0;
		simplifications += simplified.SimplifyAxB();  // replacement for simplifyLX
		// collapseXX is no longer necessary now that we're properly tracking connections

		// Handle one-connected species, notably bound solvents and metal-containing ligands.
		// TODO: check the composition of free solvents and consider connecting charged anions back to the node
		AtomSet net_1c_without_conn = simplified.GetAtoms(false).GetAtoms();
		for (AtomSet::iterator it=net_1c_without_conn.begin(); it!=net_1c_without_conn.end(); ++it) {
			ConnIntToExt unique_nbors = simplified.GetConnectedAtoms(VirtualMol(*it));
			//if (unique_nbors.size() == 1) {
			// TODO: clean up this analysis
			// Unlike the earlier algorithm, we can use the raw valence of the test point
			// because the SimplifyAxB method takes care of duplicate connections
			if ((*it)->GetValence() == 1) {
				PseudoAtom nbor_of_1c = unique_nbors.begin()->second;
				if (simplified.AtomHasRole(*it, "node")) {
					simplified.MergeAtomToAnother(*it, nbor_of_1c);
					++simplifications;
				} else if (simplified.AtomHasRole(*it, "node bridge")) {
					// probably not uncommon due to PBC and unique OBAtoms
					continue;
				} else if (simplified.AtomHasRole(*it, "linker")) {
					// Bound ligands, such as capping agents or bound solvents for ASR removal.
					// TODO: consider accounting and labeling a free solvent atom type
					// TODO: revisit this bound solvent and the earlier free solvent code
					simplified.DeleteAtomAndConns(*it);
					++simplifications;
				} else {
					obErrorLog.ThrowError(__FUNCTION__, "Unexpected atom role in the simplified net.", obWarning);
				}
			}
		}
	} while(simplifications);  // repeat until self-consistent
	//subtractMols(&mof_asr, &bound_solvent);

	// Split 4-coordinated linkers into 3+3 by convention for MIL-47, etc.
	if (mil_type_mof) {
		AtomSet for_net_4c = simplified.GetAtoms(false).GetAtoms();
		for (AtomSet::iterator it_4c=for_net_4c.begin(); it_4c!=for_net_4c.end(); ++it_4c) {
			PseudoAtom sq_4c = *it_4c;
			if (sq_4c->GetValence() == 4 && simplified.AtomHasRole(sq_4c, "linker")) {
				simplified.SplitFourVertexIntoTwoThree(sq_4c);
			}
		}
	}
	// TODO: this is where I could add other branch point detection, such as phenyl rings

	// Catenation: check that all interpenetrated nets contain identical components.
	std::vector<VirtualMol> net_components = simplified.GetAtoms().Separate();
	std::string base_formula = "";
	if (net_components.size() == 0) {
		obErrorLog.ThrowError(__FUNCTION__, "No MOFs found in the simplified net.", obError);
	} else {
		VirtualMol orig_piece = simplified.PseudoToOrig(net_components[0]);
		base_formula = orig_piece.ToOBMol(false).GetFormula();
	}
	// Compare separate topology graphs based on the orig_mol molecular formula
	for (std::vector<VirtualMol>::iterator it=net_components.begin(); it!=net_components.end(); ++it) {
		VirtualMol orig_piece = simplified.PseudoToOrig(*it);
		std::string component_formula = orig_piece.ToOBMol(false).GetFormula();
		if (component_formula != base_formula) {
			std::string err_msg =
				"Inconsistency in catenated nets.  Simplified net fragment with formula\n" +
				component_formula + " does not match first entry " + base_formula;
			obErrorLog.ThrowError(__FUNCTION__, err_msg, obWarning);
		}
	}

	// Print out the SMILES for nodes and linkers, and the detected catenation
	/*
	analysis << writeFragments(nodes.Separate(), obconv);
	analysis << writeFragments(linkers.Separate(), obconv);
	*/
	VirtualMol node_export = simplified.GetAtomsOfRole("node");
	// Handle node and node_bridge separately to match old test SMILES
	//node_export.AddVirtualMol(simplified.GetAtomsOfRole("node bridge"));
	node_export = simplified.PseudoToOrig(node_export);
	OBMol node_mol = node_export.ToOBMol();
	analysis << writeFragments(node_mol.Separate(), obconv);

	VirtualMol node_bridge_export = simplified.GetAtomsOfRole("node bridge");
	node_bridge_export = simplified.PseudoToOrig(node_bridge_export);
	OBMol node_bridge_mol = node_bridge_export.ToOBMol();
	analysis << writeFragments(node_bridge_mol.Separate(), obconv);

	VirtualMol linker_export = simplified.GetAtomsOfRole("linker");
	linker_export = simplified.PseudoToOrig(linker_export);
	OBMol linker_mol = linker_export.ToOBMol();
	analysis << writeFragments(linker_mol.Separate(), obconv);

	analysis << "Found " << net_components.size() << " simplified net(s)";


	// Write out the decomposed and simplified MOF, including bond orders
	/*
	resetBonds(&nodes);
	resetBonds(&linkers);
	writeCIF(&nodes, "Test/nodes.cif");
	writeCIF(&linkers, "Test/linkers.cif");
	*/
	writeCIF(&node_mol, "Test/nodes.cif");
	writeCIF(&linker_mol, "Test/linkers.cif");
	writeCIF(&node_bridge_mol, "Test/node_bridges.cif");

	// Write out detected solvents
	/*
	writeCIF(&free_solvent, "Test/free_solvent.cif");
	writeCIF(&bound_solvent, "Test/bound_solvent.cif");
	writeCIF(&mof_fsr, "Test/mof_fsr.cif");
	writeCIF(&mof_asr, "Test/mof_asr.cif");
	*/
	// TODO: calculate MOF-ASR and FSR separately, based on free_ and bound_solvent molecules

	// Export the simplified net
	OBMol removed_two_conn_for_topology = simplified.ToOBMol();
	writeCIF(&removed_two_conn_for_topology, "Test/removed_two_conn_for_topology.cif");
	simplified.WriteSystre("Test/topology.cgd");

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

