#include "deconstructor.h"
#include "obdetails.h"
#include "framework.h"
#include "periodic.h"
#include "topology.h"

#include <string>
#include <sstream>

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>


namespace OpenBabel
{


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



Deconstructor::Deconstructor(OBMol* orig_mof) : simplified_net(orig_mof) {
//Deconstructor::Deconstructor(OBMol* orig_mof) {
	parent_molp = orig_mof;
	// Avoid the copy constructor by using the member initializer list
	//simplified_net = Topology(parent_molp);  // TODO: does topology.h, etc., use Begin/EndModify routines?
	SetOutputDir(DEFAULT_OUTPUT_PATH);
	InitOutputFormat();
	infinite_node_detected = false;
}


void Deconstructor::InitOutputFormat() {
	// Allow overriding of obconv in derived classes
	obconv = OBConversion();
	obconv.SetOutFormat("can");  // Open Babel canonical SMILES
	obconv.AddOption("i");  // Ignore SMILES chirality for now

	/* Other potentially useful file formats, for derived classes: */
	// Universal SMILES:
	//obconv.SetOutFormat("smi");
	//obconv.AddOption("U");
	// InChI or InChIKey, with same flags as Universal SMILES:
	//obconv.SetOutFormat("inchi");
	//obconv.SetOutFormat("inchikey");
	//obconv.AddOption("F");
	//obconv.AddOption("M");
}



std::string Deconstructor::GetBasicSMILES(OBMol fragment) {
	// Get a standard, canonical SMILEs from an OBMol, e.g. to check for common fragments
	OBConversion basic_conv;
	basic_conv.SetOutFormat("can");  // Open Babel canonical SMILES
	basic_conv.AddOption("i");  // Ignore SMILES chirality for now
	return getSMILES(fragment, basic_conv);
}


void Deconstructor::DetectInitialNodesAndLinkers() {
	// Break apart the MOF to assign initial node/linker labels

	// Set up a temporary molecule for the initial metal-based fragmentation
	OBMol split_mol;
	copyMOF(parent_molp, &split_mol);

	// Find linkers by deleting bonds to metals
	// TODO: in this block, for SBU decomposition algorithms, do some manipulations to modify/restore bonds before fragmentation.
	// That will probably take the form of an optional preprocessing step before fragment assignment.
	// Similarly, there will probably be a fragmenter that breaks apart the nodes/linkers using a standard algorithm for node/linker SMILES names.
	std::vector<OBMol> fragments;
	deleteBonds(&split_mol, true);  // TODO: this could probably keep track of the deleted bond mapping
	fragments = split_mol.Separate();

	// Classify nodes and linkers based on composition.
	// Consider all single atoms and hydroxyl species as node building materials.
	std::stringstream nonmetalMsg;
	for (std::vector<OBMol>::iterator it = fragments.begin(); it != fragments.end(); ++it) {
		std::string mol_smiles = GetBasicSMILES(*it);
		VirtualMol fragment_act_atoms(parent_molp);
		fragment_act_atoms.ImportCopiedFragment(&*it);
		VirtualMol fragment_pa = simplified_net.OrigToPseudo(fragment_act_atoms);;

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
				simplified_net.DeleteAtomAndConns(*free_it, "free solvent");
			}
		} else if (it->NumAtoms() == 1) {
			nonmetalMsg << "Found a solitary atom with atomic number " << it->GetFirstAtom()->GetAtomicNum() << std::endl;
			simplified_net.SetRoleToAtoms( "node", fragment_pa);
		} else if (all_oxygens) {
			nonmetalMsg << "Found an oxygen species " << mol_smiles;
			simplified_net.SetRoleToAtoms("node", fragment_pa);  // consolidate neighboring metals in a later step
		} else {
			nonmetalMsg << "Deleting linker " << mol_smiles;
			PseudoAtom collapsed = simplified_net.CollapseFragment(fragment_pa);
			simplified_net.SetRoleToAtom("linker", collapsed);
		}
	}
	obErrorLog.ThrowError(__FUNCTION__, nonmetalMsg.str(), obDebug);
}


bool Deconstructor::CollapseNodes() {
	// Simplify all the node SBUs into single points (handling rods)
	// Returns a flag indicating whether or not an infinite, MIL-like rod was detected

	bool mil_type_mof = false;
	VirtualMol node_pa = simplified_net.GetAtomsOfRole("node");
	node_pa = simplified_net.FragmentWithIntConns(node_pa);
	std::vector<VirtualMol> node_fragments = node_pa.Separate();
	for (std::vector<VirtualMol>::iterator it=node_fragments.begin(); it!=node_fragments.end(); ++it) {
		VirtualMol fragment_mol = *it;

		OBMol fragment_obmol = fragment_mol.ToOBMol();
		if (!isPeriodicChain(&fragment_obmol)) {  // normal, nonperiodic case
			PseudoAtom collapsed = simplified_net.CollapseFragment(*it);
			simplified_net.SetRoleToAtom("node", collapsed);
		} else {  // based on sepPeriodicChains
			// TODO: consider refactoring this code to a method within topology.cpp
			// or the upcoming simplification class
			obErrorLog.ThrowError(__FUNCTION__, "Detecting infinite chains", obInfo);
			mil_type_mof = true;

			// Detect single-atom nonmetal bridging atoms
			AtomSet bridging_atoms;
			AtomSet rod_atoms = it->GetAtoms();
			for (AtomSet::iterator rod_it=rod_atoms.begin(); rod_it!=rod_atoms.end(); ++rod_it) {
				AtomSet rod_orig = simplified_net.PseudoToOrig(VirtualMol(*rod_it)).GetAtoms();
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
				simplified_net.SetRoleToAtom("node bridge", *br_it);
			}
			fragment_mol = simplified_net.FragmentWithoutConns(fragment_mol);

			if (fragment_mol.NumAtoms() == 0) {
				obErrorLog.ThrowError(__FUNCTION__, "Unexpectedly deleted all atoms in a periodic rod during simplificaiton.", obError);
				continue;
			}

			// Simplify the non-bridging metals
			std::vector<VirtualMol> rod_fragments = fragment_mol.Separate();
			for (std::vector<VirtualMol>::iterator frag_it=rod_fragments.begin(); frag_it!=rod_fragments.end(); ++frag_it) {
				if (frag_it->NumAtoms() > 1) {
					obErrorLog.ThrowError(__FUNCTION__, "Combining metal atoms within a periodic rod (likely okay, but untested code--check it).", obError);
					fragment_mol = simplified_net.FragmentWithoutConns(fragment_mol);
					PseudoAtom collapsed = simplified_net.CollapseFragment(fragment_mol);
					simplified_net.SetRoleToAtom("node", collapsed);
				}  // else, if a single-metal fragment, there's nothing to simplify
			}

		}
	}

	return mil_type_mof;
}


void Deconstructor::SimplifyTopology() {
	// Simplify the topological net (series of steps including AxB, 1-c, etc.)

	int simplifications = 0;
	do {
		simplifications = 0;
		simplifications += simplified_net.SimplifyAxB();  // replacement for simplifyLX
		// collapseXX is no longer necessary now that we're properly tracking connections

		// Handle one-connected species, notably bound solvents and metal-containing ligands.
		// TODO: check the composition of free solvents and consider connecting charged anions back to the node
		AtomSet net_1c_without_conn = simplified_net.GetAtoms(false).GetAtoms();
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

				if (simplified_net.AtomHasRole(*it, "node")) {
					if (nbor_of_1c->GetValence() == 1) {
						obErrorLog.ThrowError(__FUNCTION__, "Not collapsing 1-c node into a 1-c linker", obWarning);
					} else {
						simplified_net.MergeAtomToAnother(*it, nbor_of_1c);
						++simplifications;
					}
				} else if (simplified_net.AtomHasRole(*it, "node bridge")) {
					// probably not uncommon due to PBC and unique OBAtoms
					continue;
				} else if (simplified_net.AtomHasRole(*it, "linker")) {
					// Bound ligands, such as capping agents or bound solvents for ASR removal.
					// TODO: consider if there are cases when the bound ligand should not be removed
					simplified_net.DeleteAtomAndConns(*it, "bound solvent");
					++simplifications;
				} else {
					obErrorLog.ThrowError(__FUNCTION__, "Unexpected atom role in the simplified net.", obWarning);
				}
			}
		}
	} while(simplifications);  // repeat until self-consistent
}


int Deconstructor::CheckCatenation() {
	// Catenation: check that all interpenetrated nets contain identical components.
	// Returns the number of identified nets

	std::vector<VirtualMol> net_components = simplified_net.GetAtoms().Separate();
	std::string base_formula = "";
	if (net_components.size() == 0) {
		obErrorLog.ThrowError(__FUNCTION__, "No MOFs found in the simplified net.", obError);
	} else {
		VirtualMol orig_piece = simplified_net.PseudoToOrig(net_components[0]);
		base_formula = orig_piece.ToOBMol(false).GetFormula();
	}
	// Compare separate topology graphs based on the orig_mol molecular formula
	for (std::vector<VirtualMol>::iterator it=net_components.begin(); it!=net_components.end(); ++it) {
		VirtualMol orig_piece = simplified_net.PseudoToOrig(*it);
		std::string component_formula = orig_piece.ToOBMol(false).GetFormula();
		if (component_formula != base_formula) {
			std::string err_msg =
				"Inconsistency in catenated nets.  Simplified net fragment with formula\n" +
				component_formula + " does not match first entry " + base_formula;
			obErrorLog.ThrowError(__FUNCTION__, err_msg, obWarning);
		}
	}

	return net_components.size();
}


std::string Deconstructor::GetCatenationInfo(int num_nets) {
	// Get formatted string about catenation
	std::stringstream cat;
	cat << "Found " << num_nets << " simplified net(s)";
	return cat.str();
}


void Deconstructor::SimplifyMOF(bool write_intermediate_cifs) {
	// Runs a MOF simplfication, optionally writing intermediate CIFs

	if (write_intermediate_cifs) { WriteSimplifiedNet("simplified_test_orig.cif"); }
	DetectInitialNodesAndLinkers();
	if (write_intermediate_cifs) { WriteSimplifiedNet("test_partial.cif"); }
	infinite_node_detected = CollapseNodes();
	if (write_intermediate_cifs) { WriteSimplifiedNet("test_with_simplified_nodes.cif"); }
	SimplifyTopology();
	PostSimplification();
}


void Deconstructor::SetOutputDir(const std::string &path) {
	output_dir = path;
}


void Deconstructor::WriteCIFs() {
	// Write out accessory files: the decomposed and simplified MOF, including bond orders.
	// Also write the Systre topology file.

	WriteAtomsOfRole("node", "nodes.cif");
	WriteAtomsOfRole("linker", "linkers.cif");
	WriteAtomsOfRole("node bridge", "node_bridges.cif");

	// Write out detected solvents
	simplified_net.GetDeletedOrigAtoms("free solvent").ToCIF(GetOutputPath("free_solvent.cif"));
	simplified_net.GetDeletedOrigAtoms("bound solvent").ToCIF(GetOutputPath("bound_solvent.cif"));

	// Also write out the original MOF, desolvated into -FSR and -ASR versions
	simplified_net.PseudoToOrig(simplified_net.GetAtoms(false)).ToCIF(GetOutputPath("mof_asr.cif"));
	VirtualMol mof_fsr = simplified_net.PseudoToOrig(simplified_net.GetAtoms(false));
	mof_fsr.AddVirtualMol(simplified_net.GetDeletedOrigAtoms("bound solvent"));
	mof_fsr.ToCIF(GetOutputPath("mof_fsr.cif"));

	// Export the simplified net
	simplified_net.ToSimplifiedCIF(GetOutputPath("removed_two_conn_for_topology.cif"));
	simplified_net.WriteSystre(GetOutputPath("topology.cgd"));
}


std::string Deconstructor::GetMOFInfo() {
	// Print out the SMILES for nodes and linkers, and the detected catenation
	std::stringstream analysis;

	VirtualMol node_export = simplified_net.GetAtomsOfRole("node");
	// Handle node and node_bridge separately to match old test SMILES
	node_export = simplified_net.PseudoToOrig(node_export);
	OBMol node_mol = node_export.ToOBMol();
	analysis << writeFragments(node_mol.Separate(), obconv);

	VirtualMol node_bridge_export = simplified_net.GetAtomsOfRole("node bridge");
	node_bridge_export = simplified_net.PseudoToOrig(node_bridge_export);
	OBMol node_bridge_mol = node_bridge_export.ToOBMol();
	analysis << writeFragments(node_bridge_mol.Separate(), obconv);

	VirtualMol linker_export = simplified_net.GetAtomsOfRole("linker");
	linker_export = simplified_net.PseudoToOrig(linker_export);
	OBMol linker_mol = linker_export.ToOBMol();
	analysis << writeFragments(linker_mol.Separate(), obconv);

	analysis << GetCatenationInfo(CheckCatenation());
	return analysis.str();
}


std::string Deconstructor::GetOutputPath(const std::string &base_filename) {
	return output_dir + "/" + base_filename;
}


void Deconstructor::WriteSimplifiedNet(const std::string &base_filename) {
	simplified_net.ToSimplifiedCIF(GetOutputPath(base_filename));
}


void Deconstructor::WriteAtomsOfRole(const std::string &simplified_role, const std::string &base_filename) {
	// Writes a CIF containing the original atoms corresponding to a PA role of simplified_role.
	// If base_filename is not specified, use the simplified_role.cif

	std::string output_path = base_filename;
	if (simplified_role == "") {
		output_path = simplified_role + ".cif";
	}
	output_path = GetOutputPath(output_path);

	VirtualMol role_export = simplified_net.GetAtomsOfRole(simplified_role);
	role_export = simplified_net.PseudoToOrig(role_export);
	role_export.ToCIF(output_path);
}



MOFidDeconstructor::MOFidDeconstructor(OBMol* orig_mof) : Deconstructor(orig_mof) {
	// Note: MOFidDeconstructor would call the default constructor for Deconstructor,
	// not Deconstructor(orig_mof) unless specified above.
	// See also https://www.learncpp.com/cpp-tutorial/114-constructors-and-initialization-of-derived-classes/

	// Could initialize MOFidDeconstructor variables, etc., here
}


void MOFidDeconstructor::PostSimplification() {
	// Split 4-coordinated linkers into 3+3 by convention for MIL-47, etc.
	if (infinite_node_detected) {
		AtomSet for_net_4c = simplified_net.GetAtoms(false).GetAtoms();
		for (AtomSet::iterator it_4c=for_net_4c.begin(); it_4c!=for_net_4c.end(); ++it_4c) {
			PseudoAtom sq_4c = *it_4c;
			if (sq_4c->GetValence() == 4 && simplified_net.AtomHasRole(sq_4c, "linker")) {
				simplified_net.SplitFourVertexIntoTwoThree(sq_4c);
			}
		}
	}
	// TODO: this is where I could add other branch point detection, such as phenyl rings
}


} // end namespace OpenBabel
