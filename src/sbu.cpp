// Main MOFid code to decompose a CIF into nodes, linkers, and solvents.
// Writes the SMILES and catenation info to stdout, errors to stderr, and
// relevant CIFs and simplified topology.cgd to the Test/ directory.

// See https://openbabel.org/docs/dev/UseTheLibrary/CppExamples.html
// Get iterator help from http://openbabel.org/dev-api/group__main.shtml
// Visualize with http://baoilleach.webfactional.com/site_media/blog/emscripten/openbabel/webdepict.html

// Instead of including code to display the full MOF SMILES, just use openbabel natively:
// obabel CIFFILE -ap -ocan

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <string>
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


// Function prototypes
std::string analyzeMOF(std::string filename);
extern "C" void analyzeMOFc(const char *cifdata, char *analysis, int buflen);
extern "C" int SmilesToSVG(const char* smiles, int options, void* mbuf, unsigned int buflen);

std::string writeFragments(std::vector<OBMol> fragments, OBConversion obconv);
std::string getSMILES(OBMol fragment, OBConversion obconv);


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
	simplified.ToSimplifiedCIF("Test/simplified_test_orig.cif");

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

		if (fragment_act_atoms.GetExternalBondsOrConns().size() == 0) {
			// Assume free solvents are organic (or lone metals), so they'd be isolated without any external connections
			nonmetalMsg << "Deleting free solvent " << mol_smiles;
			AtomSet free_set = fragment_pa.GetAtoms();
			for (AtomSet::iterator free_it=free_set.begin(); free_it!=free_set.end(); ++free_it) {
				simplified.DeleteAtomAndConns(*free_it, "free solvent");
			}
		} else if (it->NumAtoms() == 1) {
			nonmetalMsg << "Found a solitary atom with atomic number " << it->GetFirstAtom()->GetAtomicNum() << std::endl;
			simplified.SetRoleToAtoms( "node", fragment_pa);
		} else if (all_oxygens) {
			nonmetalMsg << "Found an oxygen species " << mol_smiles;
			simplified.SetRoleToAtoms("node", fragment_pa);  // consolidate neighboring metals in a later step
		} else {
			nonmetalMsg << "Deleting linker " << mol_smiles;
			PseudoAtom collapsed = simplified.CollapseFragment(fragment_pa);
			simplified.SetRoleToAtom("linker", collapsed);
		}
	}
	obErrorLog.ThrowError(__FUNCTION__, nonmetalMsg.str(), obDebug);
	simplified.ToSimplifiedCIF("Test/test_partial.cif");

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

	simplified.ToSimplifiedCIF("Test/test_with_simplified_nodes.cif");

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
			// Unlike the earlier algorithm, we can use the raw valence of the test point
			// because the SimplifyAxB method takes care of duplicate connections
			if ((*it)->GetValence() == 1) {
				// Find the neighbor of the 1-coordinated PA.
				// .begin() returns the first (in this case, only) element in the internal->external map.
				PseudoAtom it_conn = VirtualMol(*it).GetExternalBondsOrConns().begin()->second;
				VirtualMol it_and_conn = VirtualMol(*it);
				it_and_conn.AddAtom(it_conn);
				PseudoAtom nbor_of_1c = it_and_conn.GetExternalBondsOrConns().begin()->second;

				if (simplified.AtomHasRole(*it, "node")) {
					simplified.MergeAtomToAnother(*it, nbor_of_1c);
					++simplifications;
				} else if (simplified.AtomHasRole(*it, "node bridge")) {
					// probably not uncommon due to PBC and unique OBAtoms
					continue;
				} else if (simplified.AtomHasRole(*it, "linker")) {
					// Bound ligands, such as capping agents or bound solvents for ASR removal.
					// TODO: consider if there are cases when the bound ligand should not be removed
					simplified.DeleteAtomAndConns(*it, "bound solvent");
					++simplifications;
				} else {
					obErrorLog.ThrowError(__FUNCTION__, "Unexpected atom role in the simplified net.", obWarning);
				}
			}
		}
	} while(simplifications);  // repeat until self-consistent

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
	VirtualMol node_export = simplified.GetAtomsOfRole("node");
	// Handle node and node_bridge separately to match old test SMILES
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
	writeCIF(&node_mol, "Test/nodes.cif");
	writeCIF(&linker_mol, "Test/linkers.cif");
	writeCIF(&node_bridge_mol, "Test/node_bridges.cif");

	// Write out detected solvents
	simplified.GetDeletedOrigAtoms("free solvent").ToCIF("Test/free_solvent.cif");
	simplified.GetDeletedOrigAtoms("bound solvent").ToCIF("Test/bound_solvent.cif");

	// Also write out the original MOF, desolvated into -FSR and -ASR versions
	simplified.PseudoToOrig(simplified.GetAtoms(false)).ToCIF("Test/mof_asr.cif");
	VirtualMol mof_fsr = simplified.PseudoToOrig(simplified.GetAtoms(false));
	mof_fsr.AddVirtualMol(simplified.GetDeletedOrigAtoms("bound solvent"));
	OBMol mof_fsr_mol = mof_fsr.ToOBMol();
	writeCIF(&mof_fsr_mol, "Test/mof_fsr.cif");

	// Export the simplified net
	simplified.ToSimplifiedCIF("Test/removed_two_conn_for_topology.cif");
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

