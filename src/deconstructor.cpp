#include "deconstructor.h"
#include "invector.h"
#include "obdetails.h"
#include "framework.h"
#include "periodic.h"
#include "topology.h"

#include <string>
#include <sstream>
#include <ostream>
#include <queue>
#include <stack>
#include <set>
#include <utility>  // std::pair

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/generic.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/elements.h>


namespace OpenBabel
{

std::set<std::string> LOGGED_ERRORS;  // global variable to keep track of reported errors in exportNormalizedMol

std::string writeFragments(std::vector<OBMol> fragments, OBConversion obconv, bool only_single_bonds) {
	// Write a list of unique SMILES for a set of fragments
	// TODO: consider stripping out extraneous tabs, etc, here or elsewhere in the code.
	std::stringstream written;
	std::set<std::string> unique_smiles;
	for (std::vector<OBMol>::iterator it = fragments.begin(); it != fragments.end(); ++it) {
		unique_smiles.insert(getSMILES(*it, obconv, only_single_bonds));  // only adds unique values in a set
	}
	for (std::set<std::string>::iterator i2 = unique_smiles.begin(); i2 != unique_smiles.end(); ++i2) {
		written << *i2;
	}

	return written.str();
}


std::string exportNormalizedMol(OBMol fragment, OBConversion obconv, bool only_single_bonds, bool unique_errors) {
	// Resets a fragment's bonding/location before format conversion
	// If only_single_bonds is set (disabled by default), only use single bonds instead of bond orders.

	// If unique_errors is set (enabled by default), only report a unique error message once per executable.
	// Otherwise, some MOFs flood the error log with warnings about aromatic bonds (raised by
	// PerceiveBondOrders within resetBonds) or unexpected valences in the InChI converter.

	std::stringstream redirected_errors;
	std::ostream* orig_err_stream = NULL;
	if (unique_errors) {
		orig_err_stream = obErrorLog.GetOutputStream();
		obErrorLog.SetOutputStream(&redirected_errors);
	}

	// A block of actual, non-error work:
	OBMol canon = fragment;
	resetBonds(&canon);
	if (only_single_bonds) {
		FOR_BONDS_OF_MOL(b, canon) {
			b->SetBondOrder(1);
		}
	}
	unwrapFragmentMol(&canon);
	std::string output_string = obconv.WriteString(&canon);

	if (unique_errors) {
		obErrorLog.SetOutputStream(orig_err_stream);  // restore the original error stream
		std::set<std::string> errors = getUniqueErrors(redirected_errors.str());
		for (std::set<std::string>::iterator it=errors.begin(); it!=errors.end(); ++it) {
			std::string err = *it;
			if (LOGGED_ERRORS.find(err) == LOGGED_ERRORS.end()) {
				LOGGED_ERRORS.insert(err);
				*orig_err_stream << err;  // re-raise the error, per the mechanism from oberror.cpp
			}
		}
	}

	return output_string;
}


std::string getSMILES(OBMol fragment, OBConversion obconv, bool only_single_bonds) {
	// Prints SMILES based on OBConversion parameters
	return exportNormalizedMol(fragment, obconv, only_single_bonds);
}


std::set<std::string> getUniqueErrors(const std::string lines_of_errors) {
	// Extract unique blocks from a string of errors

	std::set<std::string> unique_errors;
	if (lines_of_errors.length() == 0) { return unique_errors; }

	std::stringstream input_lines;
	input_lines << lines_of_errors;

	// Parse the first line and check the error format
	std::string first_line;
	std::getline(input_lines, first_line);
	if (first_line.find("====") != 0) {
		obErrorLog.ThrowError(__FUNCTION__, "Error log did not begin with ==== delimiter", obError);
		return unique_errors;
	}

	std::string curr_line;
	std::stringstream curr_block;
	curr_block << first_line << std::endl;  // initialize with the first line
	while (std::getline(input_lines, curr_line)) {
		if (curr_line.find("====") == 0) {  // new header line and block
			unique_errors.insert(curr_block.str());
			curr_block.clear();
		}
		curr_block << curr_line << std::endl;
	}
	unique_errors.insert(curr_block.str());  // finish the last block against end of string
	return unique_errors;
}


Deconstructor::Deconstructor(OBMol* orig_mof) : simplified_net(orig_mof) {
//Deconstructor::Deconstructor(OBMol* orig_mof) {
	parent_molp = orig_mof;
	points_of_extension = VirtualMol(orig_mof);
	// Avoid the Topology copy constructor by using the member initializer list
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

		if (it->NumAtoms() == 1) {
			nonmetalMsg << "Found a solitary atom with atomic number " << it->GetFirstAtom()->GetAtomicNum() << std::endl;
			simplified_net.SetRoleToAtoms( "node", fragment_pa);
		} else if (all_oxygens) {
			nonmetalMsg << "Found an oxygen species " << mol_smiles;
			simplified_net.SetRoleToAtoms("node", fragment_pa);  // consolidate neighboring metals in a later step
		} else {
			nonmetalMsg << "Deleting linker " << mol_smiles;
			simplified_net.SetRoleToAtoms("linker", fragment_pa);  // simplified later in CollapseLinkers()
		}
	}
	obErrorLog.ThrowError(__FUNCTION__, nonmetalMsg.str(), obDebug);
}


void Deconstructor::CollapseLinkers() {
	// Simplify the linker BB's into single points
	// (similar to CollapseNodes, but simpler without periodic chains)
	VirtualMol linker_pa = simplified_net.GetAtomsOfRole("linker");
	linker_pa = simplified_net.FragmentWithIntConns(linker_pa);
	std::vector<VirtualMol> linker_frags = linker_pa.Separate();
	for (std::vector<VirtualMol>::iterator it=linker_frags.begin(); it!=linker_frags.end(); ++it) {
		PseudoAtom collapsed = simplified_net.CollapseFragment(*it);
		simplified_net.SetRoleToAtom("linker", collapsed);
	}
}


bool Deconstructor::CollapseNodes() {
	// Simplify all the node SBUs into single points (handling rods)
	// Returns a flag indicating whether or not an infinite, MIL-like rod was detected

	bool mil_type_mof = false;
	VirtualMol node_pa = simplified_net.GetAtomsOfRole("node");
	node_pa = simplified_net.FragmentWithIntConns(node_pa);

	// If points of extension are defined, redefine PoE-PoE bonds as zero-atom 2-c linkers
	AtomSet full_node_pa = node_pa.GetAtoms();
	for (AtomSet::iterator it=full_node_pa.begin(); it!=full_node_pa.end(); ++it) {
		if (simplified_net.IsConnection(*it)) {
			bool connecting_poe = true;  // check if it's a PoE-PoE bond
			FOR_NBORS_OF_ATOM(nbor, *it) {
				VirtualMol nbor_pa = simplified_net.PseudoToOrig(VirtualMol(&*nbor));
				if (nbor_pa.NumAtoms() != 1) {  // PoE cannot be previously simplified
					connecting_poe = false;
				} else if (!points_of_extension.HasAtom(*(nbor_pa.GetAtoms().begin()))) {
					connecting_poe = false;
				}
			}
			// Removing the PoE-PoE connections will prevent nearby nodes from appearing as an
			// infinite chain in node_fragments (e.g. adjacent paddlewheels in hMOF-22242).
			// Since we're effectively iterating over bonds (simplified net connectors), we
			// don't have to worry about accidentally processing the same PoE-PoE bond twice.
			if (connecting_poe) {
				node_pa.RemoveAtom(*it);
				simplified_net.ConnTo2cPA(*it);
			}
		}
	}

	std::vector<VirtualMol> node_fragments = node_pa.Separate();
	for (std::vector<VirtualMol>::iterator it=node_fragments.begin(); it!=node_fragments.end(); ++it) {
		VirtualMol fragment_mol = *it;

		OBMol fragment_obmol = fragment_mol.ToOBMol();
		if (!isPeriodicChain(&fragment_obmol)) {  // normal, nonperiodic case
			PseudoAtom collapsed = simplified_net.CollapseFragment(*it);
			simplified_net.SetRoleToAtom("node", collapsed);
		} else {  // based on sepPeriodicChains
			// TODO: consider refactoring this code to a method within topology.cpp
			// or one of the deconstructor classes
			obErrorLog.ThrowError(__FUNCTION__, "Detecting infinite chains", obInfo);
			mil_type_mof = true;

			// Detect nonmetal bridging atoms (including carboxylates)
			VirtualMol bridging_atoms(fragment_mol.GetParent());
			AtomSet rod_atoms = simplified_net.FragmentWithoutConns(fragment_mol).GetAtoms();
			for (AtomSet::iterator rod_it=rod_atoms.begin(); rod_it!=rod_atoms.end(); ++rod_it) {
				AtomSet rod_orig = simplified_net.PseudoToOrig(VirtualMol(*rod_it)).GetAtoms();
				if (rod_orig.size() == 1) {
					OBAtom* single_atom = *(rod_orig.begin());
					if (!isMetal(single_atom)) {
						bridging_atoms.AddAtom(*rod_it);  // the PA, not single_atom from the original MOF
						simplified_net.SetRoleToAtom("node bridge", *rod_it);
						fragment_mol.RemoveAtom(*rod_it);
					}
				} else if (rod_orig.size() > 1) {
					obErrorLog.ThrowError(__FUNCTION__, "Unexpectedly found pre-simplified node PA with more than one original MOF atom", obError);
				}
			}

			// Reset connection PA's within the node fragment_mol
			fragment_mol = simplified_net.FragmentWithIntConns(simplified_net.FragmentWithoutConns(fragment_mol));
			if (fragment_mol.NumAtoms() == 0) {
				obErrorLog.ThrowError(__FUNCTION__, "Unexpectedly deleted all atoms in a periodic rod during simplification.", obError);
				continue;
			}

			// Simplify the non-bridging metals
			std::vector<VirtualMol> rod_fragments = fragment_mol.Separate();
			for (std::vector<VirtualMol>::iterator frag_it=rod_fragments.begin(); frag_it!=rod_fragments.end(); ++frag_it) {
				if (frag_it->NumAtoms() > 1) {
					obErrorLog.ThrowError(__FUNCTION__, "Combining metal atoms within a periodic rod (likely okay, but untested code--check it).", obError);
					PseudoAtom collapsed = simplified_net.CollapseFragment(*frag_it);
					simplified_net.SetRoleToAtom("node", collapsed);
				}  // else, if a single-metal fragment, there's nothing to simplify
			}

			// Simplify the bridging nonmetals.
			// We don't have to worry about PoE-PoE bonds in this block, because there's a
			// dedicated 2-c PA separating the PoE's in the simplified net.
			std::vector<VirtualMol> bridge_frags = simplified_net.FragmentWithIntConns(bridging_atoms).Separate();
			for (std::vector<VirtualMol>::iterator frag_it=bridge_frags.begin(); frag_it!=bridge_frags.end(); ++frag_it) {
				if (frag_it->NumAtoms() > 1) {
					obErrorLog.ThrowError(__FUNCTION__, "Combining nonmetal atoms within a periodic rod", obInfo);
					PseudoAtom collapsed = simplified_net.CollapseFragment(*frag_it);
					simplified_net.SetRoleToAtom("node bridge", collapsed);
				}  // else, if a single-atom fragment, there's nothing to simplify
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
						obErrorLog.ThrowError(__FUNCTION__, "Not collapsing 1-c node into a 1-c linker", obInfo);
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
			} else if ((*it)->GetValence() == 0) {
				// Free solvents are isolated without any external connections
				simplified_net.DeleteAtomAndConns(*it, "free solvent");
				++simplifications;
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

	if (write_intermediate_cifs) { WriteSimplifiedNet("test_simplified_orig.cif"); }
	DetectInitialNodesAndLinkers();
	CollapseLinkers();
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
	simplified_net.ToSimplifiedCIF(GetOutputPath("simplified_topology_with_two_conn.cif"));
	simplified_net.WriteSystre(GetOutputPath("topology.cgd"));
}


std::string Deconstructor::GetMOFInfo() {
	// Print out the SMILES for nodes and linkers, and the detected catenation
	std::stringstream analysis;

	const bool export_single_bonds = true;  // Using single bonds instead of bond orders for nodes

	VirtualMol node_export = simplified_net.GetAtomsOfRole("node");
	// Handle node and node_bridge separately to match old test SMILES
	node_export = simplified_net.PseudoToOrig(node_export);
	OBMol node_mol = node_export.ToOBMol();
	analysis << "# Nodes:" << std::endl;
	analysis << writeFragments(node_mol.Separate(), obconv, export_single_bonds);

	VirtualMol node_bridge_export = simplified_net.GetAtomsOfRole("node bridge");
	node_bridge_export = simplified_net.PseudoToOrig(node_bridge_export);
	OBMol node_bridge_mol = node_bridge_export.ToOBMol();
	// Considered as part of the nodes for purposes of python_smiles_parts.txt, so no subheader
	analysis << writeFragments(node_bridge_mol.Separate(), obconv, export_single_bonds);

	VirtualMol linker_export = simplified_net.GetAtomsOfRole("linker");
	linker_export = simplified_net.PseudoToOrig(linker_export);
	OBMol linker_mol = linker_export.ToOBMol();
	analysis << "# Linkers:" << std::endl;
	analysis << writeFragments(linker_mol.Separate(), obconv, !export_single_bonds);

	analysis << "# " << GetCatenationInfo(CheckCatenation());
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



MetalOxoDeconstructor::MetalOxoDeconstructor(OBMol* orig_mof) : Deconstructor(orig_mof) {
	// Note: MetalOxoDeconstructor would call the default constructor for Deconstructor,
	// not Deconstructor(orig_mof) unless specified above.
	// See also https://www.learncpp.com/cpp-tutorial/114-constructors-and-initialization-of-derived-classes/

	// Could initialize MetalOxoDeconstructor variables, etc., here
}


void MetalOxoDeconstructor::PostSimplification() {
	// Split 4-coordinated linkers into 3+3 by convention for MIL-47, etc.
	// This code was only necessary in the original MOFid deconstruction algorithm and will
	// be automatically handled in the single/all-node deconstruction algorithms.
	if (infinite_node_detected) {
		AtomSet for_net_4c = simplified_net.GetAtoms(false).GetAtoms();
		for (AtomSet::iterator it_4c=for_net_4c.begin(); it_4c!=for_net_4c.end(); ++it_4c) {
			PseudoAtom sq_4c = *it_4c;
			if (sq_4c->GetValence() == 4 && simplified_net.AtomHasRole(sq_4c, "linker")) {
				simplified_net.SplitFourVertexIntoTwoThree(sq_4c);
			}
		}
	}
}


std::vector<std::string> MetalOxoDeconstructor::PAsToUniqueInChIs(VirtualMol pa, const std::string &format) {
	// Convert PsuedoAtoms in the simplified net to their unique InChI(key) values
	// This code will strip the protonation state off of the InChIKey

	std::vector<std::string> unique_inchi;
	std::string conv_format;
	bool truncate_inchikey = false;
	if (format == "truncated inchikey") {
		conv_format = "inchikey";
		truncate_inchikey = true;
	} else if (format == "inchi" || format == "inchikey") {
		conv_format = format;
	} else {
		obErrorLog.ThrowError(__FUNCTION__, "Unexpected format flag", obError);
		return unique_inchi;
	}

	OBConversion conv;
	conv.SetOutFormat(conv_format.c_str());
	conv.AddOption("X", OBConversion::OUTOPTIONS, "SNon");  // ignoring stereochemistry, at least for now
	conv.AddOption("w");  // reduce verbosity about InChI behavior:
	// 'Omitted undefined stereo', 'Charges were rearranged', 'Proton(s) added/removed', 'Metal was disconnected'
	// See https://openbabel.org/docs/dev/FileFormats/InChI_format.html for more information.

	OBMol pa_mol = simplified_net.PseudoToOrig(pa).ToOBMol();
	std::vector<OBMol> fragments = pa_mol.Separate();
	for (std::vector<OBMol>::iterator frag=fragments.begin(); frag!=fragments.end(); ++frag) {
		std::string frag_inchi = exportNormalizedMol(*frag, conv);
		const std::string ob_newline = "\n";
		if (conv_format == "inchikey") {
			if (frag_inchi.length() != (27 + ob_newline.size())) {  // 14 + 1 + 10 + 1 + 1 + \n
				obErrorLog.ThrowError(__FUNCTION__, "Unexpected length for InChIKey", obError);
				continue;
			}
			if (truncate_inchikey) {
				// We only need the first connectivity layer if excluding stereo, isotopes, etc.
				// Per the [InChI technical FAQ](https://www.inchi-trust.org/technical-faq-2/),
				// the empty Standard InChIKey is InChIKey=MOSFIJXAXDLOML-UHFFFAOYSA-N.
				// For typical MOFs without stereo, the second layer matches "UHFFFAOYSA".
				frag_inchi = frag_inchi.substr(0, 14);
			} else {
				frag_inchi = frag_inchi.substr(0, 27);  // only stripping the final \n
			}
		}

		if (!inVector<std::string>(frag_inchi, unique_inchi)) {
			unique_inchi.push_back(frag_inchi);
		}
	}
	std::sort(unique_inchi.begin(), unique_inchi.end());  // sort alphabetically

	return unique_inchi;
}


std::string MetalOxoDeconstructor::GetMOFkey(const std::string &topology) {
	// Print out the detected MOFkey, optionally with the topology field.
	// This method is implemented in MetalOxoDeconstructor instead of the others, because the
	// organic building blocks must be intact (e.g. including carboxylates) to properly
	// calculate the MOFkey.
	std::stringstream mofkey;

	// Get unique metal atoms from the nodes
	std::vector<int> unique_elements;
	VirtualMol full_node_export = simplified_net.GetAtomsOfRole("node");
	full_node_export.AddVirtualMol(simplified_net.GetAtomsOfRole("node bridge"));
	full_node_export = simplified_net.PseudoToOrig(full_node_export);
	AtomSet node_set = full_node_export.GetAtoms();
	for (AtomSet::iterator it=node_set.begin(); it!=node_set.end(); ++it) {
		int it_element = (*it)->GetAtomicNum();
		if (isMetal(*it) && !inVector<int>(it_element, unique_elements)) {
			unique_elements.push_back(it_element);
		}
	}
	if (unique_elements.size() == 0) {
		mofkey << MOFKEY_NO_METALS;
	} else {
		std::sort(unique_elements.begin(), unique_elements.end());  // sort by atomic number
		bool first_element = true;
		for (std::vector<int>::iterator element=unique_elements.begin(); element!=unique_elements.end(); ++element) {
			if (!first_element) {
				mofkey << MOFKEY_METAL_DELIM;
			}
			mofkey << OBElements::GetSymbol(*element);
		}
	}

	// Then, write unique truncated InChIKeys (sans unnecessary layers)
	VirtualMol linker_export = simplified_net.GetAtomsOfRole("linker");
	std::vector<std::string> unique_ikeys = PAsToUniqueInChIs(linker_export, "truncated inchikey");
	if (unique_ikeys.size() == 0) {
		unique_ikeys.push_back(MOFKEY_NO_LINKERS);
	}
	for (std::vector<std::string>::iterator it=unique_ikeys.begin(); it!=unique_ikeys.end(); ++it) {
		mofkey << MOFKEY_SEP << *it;
	}

	// Add the MOFkey format signature, and topology if available
	mofkey << MOFKEY_SEP << "MOFkey-" << MOFKEY_VERSION;
	if (!topology.empty()) {
		mofkey << MOFKEY_SEP << topology;
	}

	return mofkey.str();
}


std::string MetalOxoDeconstructor::GetLinkerInChIs() {
	// Get unique, sorted InChI's for linkers in a MOF, delimited by newlines
	std::stringstream inchis;
	VirtualMol linker_export = simplified_net.GetAtomsOfRole("linker");
	std::vector<std::string> unique_inchis = PAsToUniqueInChIs(linker_export, "inchi");
	for (std::vector<std::string>::iterator it=unique_inchis.begin(); it!=unique_inchis.end(); ++it) {
		inchis << *it;  // don't need a std::endl, because Open Babel adds it by default to the export
	}
	return inchis.str();
}


std::string MetalOxoDeconstructor::GetLinkerStats(std::string sep) {
	// Get detailed stats about the linkers in a MOF, such as connectivity
	// and the SMILES-reduced InChIkey mapping

	std::stringstream output;
	std::vector<std::string> ikeys;
	std::map<std::string, std::string> ikey_to_inchi;
	std::map<std::string, std::string> ikey_to_truncated;
	std::map<std::string, std::string> ikey_to_smiles;
	std::map<std::string, std::string> ikey_to_smiles_skeleton;
	std::map<std::string, int> ikey_to_conn;
	std::map<std::string, int> ikey_to_uc_count;

	VirtualMol linker_export = simplified_net.GetAtomsOfRole("linker");
	AtomSet linker_set = linker_export.GetAtoms();
	for (AtomSet::iterator pa=linker_set.begin(); pa!=linker_set.end(); ++pa) {
		VirtualMol pa_vmol = VirtualMol(*pa);
		std::vector<std::string> pa_ikey_list = PAsToUniqueInChIs(pa_vmol, "inchikey");
		if (pa_ikey_list.size() != 1) {
			obErrorLog.ThrowError(__FUNCTION__, "Failed to calculate an InChIkey!", obWarning);
			continue;  // not updating the stats for this particular PA
		}
		std::string pa_ikey = pa_ikey_list[0];

		if (inVector<std::string>(pa_ikey, ikeys)) {
			++ikey_to_uc_count[pa_ikey];
		} else {
			ikeys.push_back(pa_ikey);
			ikey_to_uc_count[pa_ikey] = 1;
		}

		ikey_to_inchi[pa_ikey] = rtrimWhiteSpace(PAsToUniqueInChIs(pa_vmol, "inchi")[0]);
		ikey_to_truncated[pa_ikey] = rtrimWhiteSpace(PAsToUniqueInChIs(pa_vmol, "truncated inchikey")[0]);
		OBMol orig_linker = simplified_net.PseudoToOrig(pa_vmol).ToOBMol();
		const bool skeleton_flag = true;
		ikey_to_smiles[pa_ikey] = rtrimWhiteSpace(getSMILES(orig_linker, obconv, !skeleton_flag));
		ikey_to_smiles_skeleton[pa_ikey] = rtrimWhiteSpace(getSMILES(orig_linker, obconv, skeleton_flag));
		ikey_to_conn[pa_ikey] = (*pa)->GetValence();
	}

	for (std::vector<std::string>::iterator it=ikeys.begin(); it!=ikeys.end(); ++it) {
		output << *it
			<< sep << ikey_to_conn[*it]
			<< sep << ikey_to_uc_count[*it]
			<< sep << ikey_to_inchi[*it]
			<< sep << ikey_to_truncated[*it]
			<< sep << ikey_to_smiles[*it]
			<< sep << ikey_to_smiles_skeleton[*it]
			<< std::endl;
	}
	return output.str();
}



void StandardIsolatedDeconstructor::DetectInitialNodesAndLinkers() {
	// Break apart the MOF into isolated metal atoms + linkers as everything else
	// Based on the base class's Deconstructor::DetectInitialNodesAndLinkers()
	FOR_ATOMS_OF_MOL(a, *parent_molp) {
		VirtualMol a_orig(&*a);
		VirtualMol a_pseudo = simplified_net.OrigToPseudo(a_orig);
		if (isMetal(&*a)) {
			simplified_net.SetRoleToAtoms("node", a_pseudo);
		} else {
			simplified_net.SetRoleToAtoms("linker", a_pseudo);  // simplified later in CollapseLinkers()
		}
	}
}


bool StandardIsolatedDeconstructor::CollapseNodes() {
	// Nothing to do: all metal atoms are already isolated species, and we do not want to combine them
	// Return false since we're not doing anything with rod-like MOFs
	return false;
}


void StandardIsolatedDeconstructor::SimplifyTopology() {
	// Simplify the topological net adjacency matrix

	int simplifications = 0;
	do {
		simplifications = 0;

		// Check for duplicate connector sites, like the base SimplifyTopology() implementation
		simplifications += simplified_net.SimplifyAxB();

		// Simplify the adjacency matrix by outright deleting 0-c and 1-c sites
		AtomSet base_pas = simplified_net.GetAtoms(false).GetAtoms();
		for (AtomSet::iterator it=base_pas.begin(); it!=base_pas.end(); ++it) {
			if ((*it)->GetValence() == 1) {
				simplified_net.DeleteAtomAndConns(*it, "deleted 1-c site");
				++simplifications;
			} else if ((*it)->GetValence() == 0) {
				simplified_net.DeleteAtomAndConns(*it, "deleted 0-c site");
				++simplifications;
			}
		}
	} while(simplifications);  // repeat until self-consistent
}



SingleNodeDeconstructor::SingleNodeDeconstructor(OBMol* orig_mof) : Deconstructor(orig_mof) {
}


VirtualMol SingleNodeDeconstructor::GetNonmetalRingSubstituent(OBAtom* src) {
	// Traverses through neighbors (and their neighbors) using a BFS to get the
	// full set of atoms connected to src excluding rings and metals.
	// Returns a VirtualMol including src, even if neighbors are not found
	VirtualMol branch(src);
	int failed_branches = 0;  // number of branches with rings and/or metals
	std::queue<OBAtom*> to_visit;
	to_visit.push(src);

	while(!to_visit.empty()) {
		OBAtom* curr_atom = to_visit.front();
		to_visit.pop();

		FOR_NBORS_OF_ATOM(n, *curr_atom) {
			if (!branch.HasAtom(&*n) && !isMetal(&*n)) {
				if (n->IsInRing()) {
					++failed_branches;
				} else {
					to_visit.push(&*n);
				}
			}
		}
		branch.AddAtom(curr_atom);
	}

	if (failed_branches > 2) {  // allow the connections to the 2 primary ring atoms
		return VirtualMol(src);  // neighbor with specified composition does not exist
	}

	return branch;
}


std::pair<VirtualMol,VirtualMol> SingleNodeDeconstructor::CalculateNonmetalRing(OBAtom* a, OBAtom* b) {
	// Gets the nonmetal ring containing neighboring PA's a and b,
	// plus attached hydrogens and any non-ring substituents
	// (e.g. it won't grab a phenyl group or a fused ring).
	// Returns an empty molecule if a and b are not in the same ring.

	// The first VirtualMol has all ring atoms.
	// The Second VirtualMol only contains the point(s) of extension.
	// (could use a custom struct later if we need more fields)

	// Verify that a and b are neighbors
	bool found_b_nbor = false;
	FOR_NBORS_OF_ATOM(n, *a) {
		if (&*n == b) {
			found_b_nbor = true;
		}
	}
	OBMol* parent_mol = b->GetParent();
	if (!found_b_nbor) {
		obErrorLog.ThrowError(__FUNCTION__, "OBAtom's a and b must be neighbors.", obError);
		// return empty molecules
		return std::pair<VirtualMol,VirtualMol>(VirtualMol(parent_mol), VirtualMol(parent_mol));
	}

	// First run a DFS to get find the ring from a to b and get the neighbors
	const int max_ring_size = 6;
	std::stack<OBAtom*> to_visit;
	std::map<OBAtom*, OBAtom*> seen;  // <current atom, previous node>
	std::map<OBAtom*, int> seen_count;  // <current atom, size of ring>
	to_visit.push(a);
	seen[a] = NULL;
	seen_count[a] = 1;

	while (!to_visit.empty()) {
		OBAtom* curr_atom = to_visit.top();
		to_visit.pop();

		if (curr_atom == b) {
			break;  // done finding the primary a/b ring
		} else if (seen_count[curr_atom] >= max_ring_size) {
			continue;  // don't continue down nonproductive paths
		}

		FOR_NBORS_OF_ATOM(n, *curr_atom) {
			if (!isMetal(&*n) && (seen.find(&*n) == seen.end())) {
				if (curr_atom == a && &*n == b) {
					continue;  // must traverse the ring
				}
				seen[&*n] = curr_atom;
				seen_count[&*n] = seen_count[curr_atom] + 1;
				to_visit.push(&*n);
			}
		}
	}

	if (seen.find(b) == seen.end()) {
		obErrorLog.ThrowError(__FUNCTION__, "Could not find a ring containing both atoms", obError);
		return std::pair<VirtualMol,VirtualMol>(VirtualMol(parent_mol), VirtualMol(parent_mol));
	}

	// Traverse back through the ring, filling in the side substituents
	VirtualMol ring(b);
	VirtualMol atoms_with_ring_nbors(b->GetParent());  // ring atoms with nbors that are part of rings, excluding a and b
	OBAtom* curr_ring_atom = seen[b];
	while (curr_ring_atom != a) {
		ring.AddAtom(curr_ring_atom);
		if (curr_ring_atom->GetValence() > 2) {  // has non-ring coordination
			VirtualMol substituents = GetNonmetalRingSubstituent(curr_ring_atom);
			if (substituents.NumAtoms() == 1) {  // single curr_ring_atom if connection has rings, metals, etc.
				atoms_with_ring_nbors.AddAtom(curr_ring_atom);
			} else {
				ring.AddVirtualMol(substituents);
			}
		}
		curr_ring_atom = seen[curr_ring_atom];
	}
	if (atoms_with_ring_nbors.NumAtoms() > 1) {  // allow one to connect to an organic nodular BB, etc.
		obErrorLog.ThrowError(__FUNCTION__, "Possibly found a fused ring.  Skipped attachments to other rings", obWarning);
	}

	return std::pair<VirtualMol,VirtualMol>(ring, atoms_with_ring_nbors);
}


void SingleNodeDeconstructor::DetectInitialNodesAndLinkers() {
	// Identify metal-containing SBU's in the MOF
	// Build up metal SBU's based on the metals + surrounding neighbors, then assign linkers to the rest.

	// Start node identification by detecting metals
	VirtualMol nodes(parent_molp);
	FOR_ATOMS_OF_MOL(a, *parent_molp) {
		if (isMetal(&*a)) {
			nodes.AddAtom(&*a);
		}
	}

	// Add shell of nearest neighbor atoms
	AtomSet node_metals = nodes.GetAtoms();
	for (AtomSet::iterator metal=node_metals.begin(); metal!=node_metals.end(); ++metal) {
		FOR_NBORS_OF_ATOM(a, **metal) {
			nodes.AddAtom(&*a);
		}
	}

	// Process oxygen and nitrogen NN more thoroughly to determine if they should be part of the
	// node SBU or the organic linker.  Also find points of extension, like the carboxylate carbon.
	AtomSet temp_node = nodes.GetAtoms();
	VirtualMol visited_bridges(parent_molp);  // avoids counting N bridges or carboxylates twice
	for (AtomSet::iterator it=temp_node.begin(); it!=temp_node.end(); ++it) {
		PseudoAtom nn = *it;
		if (visited_bridges.HasAtom(nn)) {
			continue;  // don't visit an atom bridge twice
		}

		if (nn->GetAtomicNum() == 8) {  // oxygen
			// Check coordination environment to classify as part of the metal oxide, carboxylate, or "other"
			VirtualMol attached_hydrogens(nn->GetParent());
			OBAtom* attached_carbon = NULL;
			bool metal_oxide_or_carboxylate = true;  // only connected to metals, oxygen, and hydrogen
			FOR_NBORS_OF_ATOM(n, *nn) {
				if (!nodes.HasAtom(&*n)) {  // metals and other bound oxygens
					if (n->GetAtomicNum() == 1) {
						attached_hydrogens.AddAtom(&*n);
					} else if (n->GetAtomicNum() == 8) {
						// Like a hydrogen, the n oxygen will only get added if it's classified as a metal oxide (and not a linker)
						attached_hydrogens.AddAtom(&*n);
						obErrorLog.ThrowError(__FUNCTION__, "Found a peroxo-species.", obInfo);
						FOR_NBORS_OF_ATOM(n_nn, *n) {
							if (&*n_nn == &*nn || nodes.HasAtom(&*nn)) {
								continue;
							} else if (n_nn->GetAtomicNum() == 1) {
								attached_hydrogens.AddAtom(&*nn);
							} else if (n_nn->GetAtomicNum() == 8) {
								obErrorLog.ThrowError(__FUNCTION__, "Unexpectedly found O-O-O bond.  Classifying as linker.", obWarning);
								metal_oxide_or_carboxylate = false;
							} else {
								obErrorLog.ThrowError(__FUNCTION__, "Unexpectedly found O-O-X bond.  Classifying as linker.", obWarning);
								metal_oxide_or_carboxylate = false;
							}
						}
					} else if (n->GetAtomicNum() == 6) {
						if (attached_carbon) {
							obErrorLog.ThrowError(__FUNCTION__, "Found a metal-bound oxygen with multiple carbon neighbors!  Results may have inconsistencies.", obError);
						}
						attached_carbon = &*n;
					} else {  // oxygen nn is bonded to another nonmetal n
						metal_oxide_or_carboxylate = false;
					}
				}
			}

			if (!metal_oxide_or_carboxylate) {
				// Oxygen as part of a non-carboxylate, nonmetal cluster
				nodes.RemoveAtom(nn);  // not part of an SBU cluster
			} else if (!attached_carbon) {  // metal oxide
				// Keep the oxygen and its bound H's (or peroxide)
				nodes.AddVirtualMol(attached_hydrogens);
			} else {  // carboxylate
				// Check the valence of the carbon for carboxylate chemistry
				bool is_carboxylate = true;
				int carbon_c_nbors = 0;
				int carbon_o_nbors = 0;
				OBAtom* other_o_nbor = NULL;
				FOR_NBORS_OF_ATOM(c_nn, attached_carbon) {
					if (c_nn->GetAtomicNum() == 6) {
						++carbon_c_nbors;
					} else if (c_nn->GetAtomicNum() == 8) {
						++carbon_o_nbors;
						if (&*c_nn != nn) {
							other_o_nbor = &*c_nn;
						}
					} else {
						is_carboxylate = false;
					}
				}
				if (carbon_c_nbors != 1 || carbon_o_nbors != 2) {
					is_carboxylate = false;
				}

				if (!is_carboxylate) {
					nodes.RemoveAtom(nn);  // nn oxygen is not actually part of a carboxylate
				} else {
					points_of_extension.AddAtom(attached_carbon);
					VirtualMol carboxylate(nn);
					carboxylate.AddAtom(attached_carbon);
					carboxylate.AddAtom(other_o_nbor);
					visited_bridges.AddVirtualMol(carboxylate);
					nodes.AddVirtualMol(carboxylate);

					// Get hydrogens attached to the second oxygen atom
					FOR_NBORS_OF_ATOM(b_h, *other_o_nbor) {
						if (b_h->GetAtomicNum() == 1) {
							attached_hydrogens.AddAtom(&*b_h);
						}
					}
					// Bound H's share the same fate/assignment as their oxygen atoms
					nodes.AddVirtualMol(attached_hydrogens);
				}
			}

		} else if (nn->GetAtomicNum() == 7) {  // nitrogen
			// Check for M-N-N-M (if the nitrogen a neighbor that's a NN to a metal)
			OBAtom* bridge_bw_metals = NULL;
			FOR_NBORS_OF_ATOM(n, *nn) {
				if (n->GetAtomicNum() == 7 && nodes.HasAtom(&*n)) {
					// FIXME: we should more generally handle cases other than N-N binding, without breaking ZIFs.
					// Probably implement by checking bridge_bw_metals only for a missing NN, then maybe
					// constructing a new test molecule via the nn_ring and existing nodes for periodicity?
					// TODO: also considering cases like hMOF-17399 so we don't grab the nitrogen as part of the node?
					// It's increasingly sounding like we should check for periodicity as well as number of connections to the same metal cluster
					bridge_bw_metals = &*n;
					// There is no problem with an N-N-N ring, e.g. MFU-4L, so deleting this code:
					//if (bridge_bw_metals) {  // already defined
					//	obErrorLog.ThrowError(__FUNCTION__, "Unexpectedly found N-N-N species in node.", obWarning);
					//}
				}
			}
			if (!bridge_bw_metals) {  // Binding end-on, e.g. pillared MOFs, porphyrins, and probably amines
				nodes.RemoveAtom(nn);  // not part of an SBU cluster
			} else {  // If there is a bridge, grab the whole ring
				std::pair<VirtualMol,VirtualMol> ring_info = CalculateNonmetalRing(nn, bridge_bw_metals);
				VirtualMol nn_ring = ring_info.first;
				nodes.AddVirtualMol(nn_ring);
				visited_bridges.AddVirtualMol(nn_ring);

				VirtualMol nn_poe = ring_info.second;
				// But don't classify ring atoms bound to the metals as PoE's
				AtomSet nn_poe_set = nn_poe.GetAtoms();
				for (AtomSet::iterator poe_it=nn_poe_set.begin(); poe_it!=nn_poe_set.end(); ++poe_it) {
					if (nodes.HasAtom(*poe_it)) {  // FIXME: not exactly correct: we've already added the nodes!
						// TODO: We could consider implementing via a neighbor search, but it still doesn't
						// account for the PoE-PoE deletion, etc., for fused rings.
						// Or, maybe we just accept that PoE will not be correctly identify if rings are fused.
						// In that case, it's not a true PoE anyway
						nn_poe.RemoveAtom(*poe_it);
					}
				}
				points_of_extension.AddVirtualMol(nn_poe);
			}
		}
	}

	// Assign node PA's
	VirtualMol node_pa = simplified_net.OrigToPseudo(nodes);
	simplified_net.SetRoleToAtoms("node", node_pa);

	// Assign linker type to PA's
	// The remainder of the atoms are either organic SBUs or 2-c linkers by definition.
	// Keep points of extension as part of the linkers to avoid issues with PoE-PoE bonds.
	// Incorporating them as linkers will automatically take care of details like 2-c sites.
	VirtualMol orig_linkers(parent_molp);
	FOR_ATOMS_OF_MOL(a, *parent_molp) {
		if (!nodes.HasAtom(&*a)) {
			orig_linkers.AddAtom(&*a);
		}
	}
	VirtualMol pa_linkers = simplified_net.OrigToPseudo(orig_linkers);
	simplified_net.SetRoleToAtoms("linker", pa_linkers);

	// TODO: consider rewriting the base MOFid separation algorithm to be additive like this one
}


void SingleNodeDeconstructor::WriteCIFs() {
	// Call base class exporter, plus the new outputs
	Deconstructor::WriteCIFs();
	points_of_extension.ToCIF(GetOutputPath("points_of_extension.cif"));
	WriteSBUs("node_sbus_no_ext_conns.cif", false, false);
	WriteSBUs("node_sbus_with_ext_conns.cif", true, true);
}


void SingleNodeDeconstructor::WriteSBUs(const std::string &base_filename, bool external_bond_pa, bool external_conn_pa) {
	// Write out the combined SBU's, including points of extension.
	// If external_bond_pa is true, also write out pseudoatoms corresponding to the next
	// external atom connected to the point of extension.
	// If external_conn_pa is also true, also write out psuedoatoms for other connections
	// from the SBU, not just the point of extension.

	// Copy PoE's as new atoms
	VirtualMol sbus = simplified_net.PseudoToOrig(simplified_net.GetAtomsOfRole("node"));
	sbus.AddVirtualMol(simplified_net.PseudoToOrig(simplified_net.GetAtomsOfRole("node bridge")));
	MappedMol sbu_mapping;
	sbus.CopyToMappedMol(&sbu_mapping);
	OBMol* sbu_molp = &(sbu_mapping.mol_copy);

	// Find PoE-PoE bonds then delete them
	FOR_BONDS_OF_MOL(b, *parent_molp) {
		OBAtom* begin = b->GetBeginAtom();
		OBAtom* end = b->GetEndAtom();
		if (points_of_extension.HasAtom(begin) && points_of_extension.HasAtom(end)) {
			OBAtom* mol_begin = sbu_mapping.origin_to_copy[begin];
			OBAtom* mol_end = sbu_mapping.origin_to_copy[end];
			sbu_molp->DeleteBond(sbu_molp->GetBond(mol_begin, mol_end));
		}
	}

	// Form external bonds to POE's and otherwise
	ConnIntToExt ext_bonds_orig;  // mapping of bonds to form, based on the original MOF OBMol
	std::map<AtomPair, int> ext_bond_elements;  // which element (atomic number) to use for the PA

	if (external_bond_pa) {
		AtomSet poe_set = points_of_extension.GetAtoms();
		for (AtomSet::iterator it=poe_set.begin(); it!=poe_set.end(); ++it) {
			OBAtom* poe_orig_atom = *it;  // PoE in the original MOF OBMol
			FOR_NBORS_OF_ATOM(nbor, *poe_orig_atom) {
				// Check PoE-external bonds or PoE-PoE connections.
				// It's okay to have both PoE1-PoE2 and PoE2-PoE1 pairs, because
				// we want both directions to be represented in the external conns.
				if (!sbus.HasAtom(&*nbor) || points_of_extension.HasAtom(&*nbor)) {
					AtomPair bond_pair(poe_orig_atom, &*nbor);
					ext_bonds_orig.insert(bond_pair);
					ext_bond_elements[bond_pair] = POE_EXTERNAL_ELEMENT;
				}
			}
		}
	}
	if (external_conn_pa) {
		// Similar to external_bond_pa, but bonding over the SBU
		AtomSet sbu_set = sbus.GetAtoms();
		for (AtomSet::iterator it=sbu_set.begin(); it!=sbu_set.end(); ++it) {
			OBAtom* sbu_orig_atom = *it;  // SBU atom in the original MOF OBMol
			FOR_NBORS_OF_ATOM(nbor, *sbu_orig_atom) {
				if (!sbus.HasAtom(&*nbor)) {
					AtomPair bond_pair(sbu_orig_atom, &*nbor);
					// Not overwiting existing pairs from external_bond_pa
					if (ext_bonds_orig.find(bond_pair) == ext_bonds_orig.end()) {
						ext_bonds_orig.insert(bond_pair);
						ext_bond_elements[bond_pair] = SBU_EXTERNAL_ELEMENT;
					}
				}
			}
		}
	}

	// Connect the PoE or other internal SBU atom to the external
	OBUnitCell* uc = getPeriodicLattice(parent_molp);
	for (ConnIntToExt::iterator it=ext_bonds_orig.begin(); it!=ext_bonds_orig.end(); ++it) {
		int connection_element = ext_bond_elements[*it];
		PseudoAtom int_orig_atom = it->first;
		PseudoAtom ext_orig_atom = it->second;

		// Calculate a distance 1/3 of the way from the PoE/internal to the external atom to avoid
		// overlapping atoms when two POEs are next to each other in the MOF (e.g. POXHER_clean.cif)
		vector3 int_pos = int_orig_atom->GetVector();
		vector3 ext_pos = uc->UnwrapCartesianNear(ext_orig_atom->GetVector(), int_pos);
		vector3 conn_pos = uc->WrapCartesianCoordinate((2.0*int_pos + ext_pos) / 3.0);

		PseudoAtom new_conn_atom = formAtom(sbu_molp, conn_pos, connection_element);
		OBAtom* int_mol_atom = sbu_mapping.origin_to_copy[int_orig_atom];
		formBond(sbu_molp, int_mol_atom, new_conn_atom, 1);
	}

	writeCIF(sbu_molp, GetOutputPath(base_filename));
}



AllNodeDeconstructor::AllNodeDeconstructor(OBMol* orig_mof) : SingleNodeDeconstructor(orig_mof) {
	branches_orig = VirtualMol(orig_mof);
	branch_points_orig = VirtualMol(orig_mof);

	OBMol* simplified_molp = simplified_net.GetAtoms().GetParent();  // temporary var
	branches_pa = initMOFwithUC(simplified_molp);
	branch_points_pa = initMOFwithUC(simplified_molp);
}

void AllNodeDeconstructor::CollapseLinkers() {
	// Detect branch points in the linkers, then collapse into PA's
	// Overrides a base class version, e.g. SingleNodeDeconstructor::CollapseLinkers() or inherited

	// Start with a description from the single-node decomposition
	VirtualMol linker_pa = simplified_net.GetAtomsOfRole("linker");
	linker_pa = simplified_net.FragmentWithIntConns(linker_pa);
	std::vector<VirtualMol> single_node_frags = linker_pa.Separate();

	for (std::vector<VirtualMol>::iterator frag=single_node_frags.begin(); frag!=single_node_frags.end(); ++frag) {
		MappedMol frag_all_node;  // simplified all-node fragment
		frag->CopyToMappedMol(&frag_all_node);

		// Get connection sites from fragment to external PA's
		ConnIntToExt pa_int_to_ext = frag->GetExternalBondsOrConns();  // already has IntConns from linker_pa assignment
		VirtualMol int_conns(frag->GetParent());
		for (ConnIntToExt::iterator it=pa_int_to_ext.begin(); it!=pa_int_to_ext.end(); ++it) {
			int_conns.AddAtom(it->first);
		}

		// Identify branches and branch points
		TreeDecomposition(&frag_all_node, int_conns);

		// Collapse linkers by applying the changes indicated from frag_all_node.
		// Branch points are kept, and branches (internal/connectors) can are simplified into bonds.
		// Recall that collapseSBU will automatically handle external bonds to the rest of the simplified net.

		// One subtlety is that we needed the ConnectionTable PA's to keep everything glued together for TreeDecomposition.
		// However, when we start modifying/simplifying the simplified net, the table of these connectors will change,
		// causing inconsistencies in the copy_pa_to_multiple mapping (former connection links in copy_pa_to_multiple will
		// reference PA's that no longer exist in the simplified net, as bonds are formed and broken).
		// To avoid those issues, remove the connectors now.  CollapseFragment and related methods in the Topology class
		// are smart enough to grab the correct, fully updated pointers for the ConnectionTable class.
		// The case of a connector-only PA will be handled in the branch simplification code below.
		for (std::map<PseudoAtom,VirtualMol>::iterator it=frag_all_node.copy_pa_to_multiple.begin(); it!=frag_all_node.copy_pa_to_multiple.end(); ++it) {
			if (it->second.NumAtoms() == 0) {
				obErrorLog.ThrowError(__FUNCTION__, "Found empty fragment before connection removal", obWarning);
			}
			it->second = simplified_net.FragmentWithoutConns(it->second);
		}

		// Collapse branch points, keeping track of the new PA's for branch simplification later
		std::map<PseudoAtom, PseudoAtom> mapped_to_net;
		std::map<PseudoAtom, PseudoAtom> net_to_mapped;
		FOR_ATOMS_OF_MOL(pa, frag_all_node.mol_copy) {
			if (pa->GetAtomicNum() == TREE_BRANCH_POINT) {
				VirtualMol pa_network = frag_all_node.copy_pa_to_multiple[&*pa];
				PseudoAtom collapsed = simplified_net.CollapseFragment(pa_network);
				simplified_net.SetRoleToAtom("linker", collapsed);
				branch_points_orig.AddVirtualMol(simplified_net.PseudoToOrig(collapsed));
				formAtom(&branch_points_pa, collapsed->GetVector(), TREE_BRANCH_POINT);

				mapped_to_net[&*pa] = collapsed;
				net_to_mapped[collapsed] = &*pa;
			} else if (pa->GetAtomicNum() != TREE_EXT_CONN && pa->GetAtomicNum() != TREE_INT_BRANCH) {
				obErrorLog.ThrowError(__FUNCTION__, "Found unexpected pseudoatom type in MappedMol", obError);
			}
		}

		// Collapse branches as bonds
		FOR_ATOMS_OF_MOL(pa, frag_all_node.mol_copy) {
			if (pa->GetAtomicNum() == TREE_EXT_CONN || pa->GetAtomicNum() == TREE_INT_BRANCH) {
				// Check expected valencies
				if (pa->GetAtomicNum() == TREE_EXT_CONN) {
					if (pa->GetValence() != 1) {
						obErrorLog.ThrowError(__FUNCTION__, "Found TREE_EXT_CONN with an unexpected implicit valence", obError);
						continue;
					}
				} else if (pa->GetAtomicNum() == TREE_INT_BRANCH) {
					if (pa->GetValence() != 2) {
						obErrorLog.ThrowError(__FUNCTION__, "Found TREE_INT_BRANCH with an unexpected valence", obError);
						continue;
					}
				}

				std::vector<PseudoAtom> real_nbors_from_map;
				FOR_NBORS_OF_ATOM(nbor, *pa) {
					if (nbor->GetAtomicNum() == TREE_BRANCH_POINT) {
						real_nbors_from_map.push_back(mapped_to_net[&*nbor]);
					} else {
						obErrorLog.ThrowError(__FUNCTION__, "Found a branch with a non-branch-point neighbor.  Skipping branch simplification", obError);
						continue;
					}
				}

				// Collapse PA's in the simplified net (unless it's a ghost branch from a connector PA)
				VirtualMol pa_network = frag_all_node.copy_pa_to_multiple[&*pa];
				if (pa_network.NumAtoms() == 0) {
					obErrorLog.ThrowError(__FUNCTION__, "Found empty branch (it was formerly a connector PA, unless this message is preceded by an empty fragment warning from above).", obInfo);
					if (pa->GetAtomicNum() == TREE_EXT_CONN) {
						obErrorLog.ThrowError(__FUNCTION__, "Unexpectedly found a TREE_EXT_CONN mapped to zero atoms", obError);
					}
					// No work is needed, since the simplified net already has a zero-atom connector PA between the two branch points (without intermediate real PA's)
					continue;
				}
				PseudoAtom collapsed = simplified_net.CollapseFragment(pa_network);
				simplified_net.SetRoleToAtom("linker", collapsed);
				branches_orig.AddVirtualMol(simplified_net.PseudoToOrig(collapsed));
				formAtom(&branches_pa, collapsed->GetVector(), pa->GetAtomicNum());

				std::map<PseudoAtom, vector3> external_nbors_to_restore;  // <PA in simplified net, conn location>
				VirtualMol external_conns_to_delete(frag->GetParent());
				FOR_NBORS_OF_ATOM(nbor, *collapsed) {
					PseudoAtom curr_conn = &*nbor;
					PseudoAtom curr_nbor = simplified_net.GetOtherEndpoint(curr_conn, collapsed);
					if (!inVector<PseudoAtom>(curr_nbor, real_nbors_from_map)) {
						external_nbors_to_restore[curr_nbor] = curr_conn->GetVector();
						external_conns_to_delete.AddAtom(curr_conn);
					}
				}
				// If it's an internal branch, we're only merging to one of the two ends, so be sure to handle the second end
				if (pa->GetAtomicNum() == TREE_INT_BRANCH) {
					PseudoAtom other_end = real_nbors_from_map[1];
					PseudoAtom other_conn = NULL;

					FOR_NBORS_OF_ATOM(n_conn, *collapsed) {
						if (simplified_net.GetOtherEndpoint(&*n_conn, collapsed) == other_end) {
							if (other_conn) {
								obErrorLog.ThrowError(__FUNCTION__, "Found multiple matching connections for real_nbors_from_map[1].", obWarning);
							}
							other_conn = &*n_conn;
						}
					}
					if (!other_conn) {
						obErrorLog.ThrowError(__FUNCTION__, "Could not find a matching connection for real_nbors_from_map[1].", obError);
						continue;
					}
					external_nbors_to_restore[other_end] = other_conn->GetVector();
					external_conns_to_delete.AddAtom(other_conn);
				}
				AtomSet external_conns_set = external_conns_to_delete.GetAtoms();
				for (AtomSet::iterator it=external_conns_set.begin(); it!=external_conns_set.end(); ++it) {
					simplified_net.DeleteConnection(*it);
				}
				// Make the connection to the only mapped nbor (TREE_EXT_CONN) or
				// randomly one of the two mapped nbors (TREE_INT_BRANCH)
				simplified_net.MergeAtomToAnother(collapsed, real_nbors_from_map[0]);

				// Restore external bonds formerly attached to the branch
				for (std::map<PseudoAtom,vector3>::iterator it=external_nbors_to_restore.begin(); it!=external_nbors_to_restore.end(); ++it) {
					simplified_net.ConnectAtoms(real_nbors_from_map[0], it->first, &(it->second));
				}
			}
		}

	}
}


void AllNodeDeconstructor::TreeDecomposition(MappedMol *fragment_to_simplify, VirtualMol connection_points) {
	// Collapse an abstract, organic fragment molecule into its "all-node" representation.
	// Inspired by Algorithm 2 of arXiv:1802.04364v2, except instead of calculating the minimum
	// spanning tree, our goal is to fully simplify the graph and account for endpoint connectors.

	// Keeps connection_points from fragment_to_simplify->origin_molp (the simplified topology OBMol)
	// as branches instead of collapsing them into branch points.  In the end, there will be 3 atom types:
	// 1. TREE_BRANCH_POINT: internal branch point, i.e. the key difference between the all-node and single-node algorithms
	// 2. TREE_EXT_CONN: branch referencing connection_points, which link externally to the rest of the framework
	// 3. TREE_INT_BRANCH: internal branches connecting two branch points

	// My apologies for the long function:

	// Option to write out the intermediate results as CIFs for debugging or demonstrating the algorithm.
	// Note: TreeDecomposition is usually called many times per MOF, so only the last set of CIFs will be kept
	const bool DEBUG_WITH_CIFS = false;

	// Aliases
	MappedMol* frag_map = fragment_to_simplify;
	OBMol* frag_molp = &(frag_map->mol_copy);

	// Invalidate 1:1 mappings since we will be collapsing PA's
	frag_map->origin_to_copy.clear();
	frag_map->copy_to_origin.clear();


	// If there are no external connection points (e.g. free solvent), simplify to a single point
	if (connection_points.NumAtoms() == 0) {
		vector3 centroid = getCentroid(frag_molp, false);
		VirtualMol frag_atoms(frag_molp);
		FOR_ATOMS_OF_MOL(a, *frag_molp) {
			frag_atoms.AddAtom(&*a);
		}
		AtomSet frag_atoms_set = frag_atoms.GetAtoms();

		PseudoAtom single_point = formAtom(frag_molp, centroid, TREE_BRANCH_POINT);
		VirtualMol all_origins(frag_map->origin_molp);
		for (AtomSet::iterator it=frag_atoms_set.begin(); it!=frag_atoms_set.end(); ++it) {
			all_origins.AddVirtualMol(frag_map->copy_pa_to_multiple[*it]);
			frag_map->copy_pa_to_multiple.erase(*it);
			frag_molp->DeleteAtom(*it);
		}
		frag_map->copy_pa_to_multiple[single_point] = all_origins;

		if (frag_map->copy_pa_to_multiple.size() != 1) {
			obErrorLog.ThrowError(__FUNCTION__, "Assertion: incorrect number of entries in MappedMol.copy_pa_to_multiple", obError);
		}
		return;
	}


	// Rings to collapsed pseudo atoms
	if (DEBUG_WITH_CIFS) {writeCIF(frag_molp, GetOutputPath("debug_tree_1_initial.cif")); }
	CollapseRings(frag_map);
	if (DEBUG_WITH_CIFS) {writeCIF(frag_molp, GetOutputPath("debug_tree_2_collapsed_rings.cif")); }


	// Label connections to external atoms, including rings which could also be internal branch points.
	// The valence of connection sites will be reevaluated later, considering external connection points as
	// additional valence in the graph vertices.
	// (e.g. differentiating "2-c" ring in bispyridine vs. a ">=3-c" ring branch point with internal bonds as well)
	FOR_ATOMS_OF_MOL(pa, *frag_molp) {
		AtomSet orig_from_pa = frag_map->copy_pa_to_multiple[&*pa].GetAtoms();
		if (orig_from_pa.size() == 0) {
			obErrorLog.ThrowError(__FUNCTION__, "Unexpectedly found a PA without origin atoms during connection labeling", obError);
		}
		int num_conn_sites = 0;
		for (AtomSet::iterator it=orig_from_pa.begin(); it!=orig_from_pa.end(); ++it) {
			if (connection_points.HasAtom(*it)) { ++num_conn_sites; }
		}
		if (num_conn_sites == 1) {
			pa->SetAtomicNum(TREE_EXT_CONN);
		} else if (num_conn_sites > 1) {
			// Technically, this rule is overly broad and will also classify simple 2-c linkers as
			// BP's (even if they don't contain an additional internal connection), but it's
			// otherwise a simple method to detect branch point rings with 2+ external connections
			// which may also include an internal bond.
			pa->SetAtomicNum(TREE_BRANCH_POINT);
		}
	}
	if (DEBUG_WITH_CIFS) {writeCIF(frag_molp, GetOutputPath("debug_tree_3_label_conns.cif")); }


	// Start simplifying the fragment
	// First, iteratively collapse 1-c until self-consistent
	int simplifications = 0;
	do {
		simplifications = 0;
		ConnIntToExt pa_to_1c;
		FOR_ATOMS_OF_MOL(a, *frag_molp) {
			if (a->GetValence() == 1 && a->GetAtomicNum() != TREE_EXT_CONN && a->GetAtomicNum() != TREE_BRANCH_POINT) {
				FOR_NBORS_OF_ATOM(nbor, *a) {  // get the 1 neighbor (inner to fragment)
					pa_to_1c.insert(AtomPair(&*nbor, &*a));
				}
			}
		}

		// Apply changes to the fragment
		for (ConnIntToExt::iterator it=pa_to_1c.begin(); it!=pa_to_1c.end(); ++it) {
			PseudoAtom pa_kept = it->first;
			PseudoAtom pa_1c = it->second;

			VirtualMol kept_vmol = frag_map->copy_pa_to_multiple[pa_kept];
			kept_vmol.AddVirtualMol(frag_map->copy_pa_to_multiple[pa_1c]);
			frag_map->copy_pa_to_multiple[pa_kept] = kept_vmol;
			frag_map->copy_pa_to_multiple.erase(pa_1c);

			frag_molp->DeleteAtom(pa_1c);
			++simplifications;
		}
	} while (simplifications);
	if (DEBUG_WITH_CIFS) {writeCIF(frag_molp, GetOutputPath("debug_tree_4_simplify_1c.cif")); }


	// At this point, the valence of external connectors and internal branch points will be set.
	// Differentiate between truly external connectors and internal branch points containing connection sites.
	FOR_ATOMS_OF_MOL(a, frag_molp) {
		if (a->GetAtomicNum() == TREE_EXT_CONN) {
			if (a->GetValence() == 0) {
				if (frag_molp->NumAtoms() > 1) {
					obErrorLog.ThrowError(__FUNCTION__, "AssertionError: found TREE_EXT_CONN PA disconnected from the rest of the graph", obError);
				}
				a->SetAtomicNum(TREE_BRANCH_POINT);  // everything collapsed to a single point, e.g. pyrazine
			} else if (a->GetValence() == 1) {
				// Keep it as type TREE_EXT_CONN: 2-c with explicit connection plus an implicit external connection
			} else {  // total valence is 3+
				a->SetAtomicNum(TREE_BRANCH_POINT);
			}
		}
	}
	if (DEBUG_WITH_CIFS) {writeCIF(frag_molp, GetOutputPath("debug_tree_5_convert_conns_to_bp.cif")); }


	// Collapse external branches into the PA's for external connectors
	// TODO: consider refactoring this loop and the previous one into a more general subroutine?
	// WARNING: this loop simplifies in the opposite order of the previous loop: PA collapses into the connector
	do {
		simplifications = 0;
		ConnIntToExt pa_to_conn;
		FOR_ATOMS_OF_MOL(a, *frag_molp) {
			if (a->GetAtomicNum() == TREE_EXT_CONN) {
				FOR_NBORS_OF_ATOM(nbor, *a) {
					// nbor should only connects to TREE_EXT_CONN and (optionally) one other site,
					// otherwise it's a branch point.
					// Handle the conn-conn case separately as another step
					if (nbor->GetValence() < 3 && nbor->GetAtomicNum() != TREE_EXT_CONN && nbor->GetAtomicNum() != TREE_BRANCH_POINT) {
						// Only allow one connector to claim the nbor.  Otherwise, you could have
						// two connectors competing to simplify conn-2c-conn.
						bool new_pa_simplification = true;
						for (ConnIntToExt::iterator it=pa_to_conn.begin(); it!=pa_to_conn.end(); ++it) {
							if (it->first == &*nbor) {
								new_pa_simplification = false;
							}
						}
						if (new_pa_simplification) {
							pa_to_conn.insert(AtomPair(&*nbor, &*a));
						}
					}
				}
			}
		}

		// Apply changes to the fragment
		for (ConnIntToExt::iterator it=pa_to_conn.begin(); it!=pa_to_conn.end(); ++it) {
			PseudoAtom pa_simplified = it->first;
			PseudoAtom pa_conn = it->second;

			VirtualMol conn_vmol = frag_map->copy_pa_to_multiple[pa_conn];
			conn_vmol.AddVirtualMol(frag_map->copy_pa_to_multiple[pa_simplified]);
			frag_map->copy_pa_to_multiple[pa_conn] = conn_vmol;
			frag_map->copy_pa_to_multiple.erase(pa_simplified);

			FOR_NBORS_OF_ATOM(nbor, pa_simplified) {
				if (&*nbor != pa_conn) {
					formBond(frag_molp, pa_conn, &*nbor, 1);
				}
			}

			frag_molp->DeleteAtom(pa_simplified);
			++simplifications;
		}
	} while (simplifications);
	if (DEBUG_WITH_CIFS) {writeCIF(frag_molp, GetOutputPath("debug_tree_6_collapse_conns.cif")); }


	// Combine adjacent connection sites having one explicit bond
	// (e.g. the two carboxylates of BDC in MOF-5) into branch points
	std::vector<VirtualMol> pairs_of_1c_conns;  // can't use a set because VirtualMols have no comparison operator
	VirtualMol simplified_1c(frag_molp);
	FOR_ATOMS_OF_MOL(a, frag_molp) {
		if (a->GetAtomicNum() == TREE_EXT_CONN && a->GetValence() == 1) {
			FOR_NBORS_OF_ATOM(a_nbor, *a) {
				if (a_nbor->GetAtomicNum() == TREE_EXT_CONN && a_nbor->GetValence() == 1) {
					if (!simplified_1c.HasAtom(&*a) && !simplified_1c.HasAtom(&*a_nbor)) {
						simplified_1c.AddAtom(&*a);
						simplified_1c.AddAtom(&*a_nbor);
						VirtualMol pair_conns(frag_molp);
						pair_conns.AddAtom(&*a);
						pair_conns.AddAtom(&*a_nbor);
						pairs_of_1c_conns.push_back(pair_conns);
					}
				}
			}
		}
	}
	for (std::vector<VirtualMol>::iterator it=pairs_of_1c_conns.begin(); it!=pairs_of_1c_conns.end(); ++it) {
		// Make a new PA at the midpoint of the connection sites
		VirtualMol pair_vmol = *it;
		OBMol pair_mol = pair_vmol.ToOBMol();
		vector3 pair_centroid = getCentroid(&pair_mol, false);
		PseudoAtom new_bp = formAtom(frag_molp, pair_centroid, TREE_BRANCH_POINT);

		// Update accounting in MappedMol.  No bonds to reform, because both PA's were 1c to each other
		VirtualMol new_origin_map = VirtualMol(frag_map->origin_molp);
		AtomSet pair_atoms = pair_vmol.GetAtoms();
		for (AtomSet::iterator pair_atom=pair_atoms.begin(); pair_atom!=pair_atoms.end(); ++pair_atom) {
			new_origin_map.AddVirtualMol(frag_map->copy_pa_to_multiple[*pair_atom]);
			frag_map->copy_pa_to_multiple.erase(*pair_atom);
			frag_molp->DeleteAtom(*pair_atom);
		}
		frag_map->copy_pa_to_multiple[new_bp] = new_origin_map;
	}
	if (DEBUG_WITH_CIFS) {writeCIF(frag_molp, GetOutputPath("debug_tree_7_adjacent_conns.cif")); }


	// Simplify remaining 2-c sites, which are internal branches (e.g. alkyne spacers)
	// WARNING: copy-pasted with some modifications from previous block
	std::vector<VirtualMol> branches_to_combine;
	VirtualMol visited_branches(frag_molp);
	FOR_ATOMS_OF_MOL(a, frag_molp) {
		if (a->GetValence() == 2 && a->GetAtomicNum() != TREE_EXT_CONN && a->GetAtomicNum() != TREE_BRANCH_POINT && !visited_branches.HasAtom(&*a)) {
			// Use a BFS to get all adjacent 2-c sites
			VirtualMol a_nbors(frag_molp);
			std::queue<PseudoAtom> to_visit;
			to_visit.push(&*a);

			while (!to_visit.empty()) {
				PseudoAtom curr_bfs = to_visit.front();
				to_visit.pop();
				visited_branches.AddAtom(curr_bfs);
				a_nbors.AddAtom(curr_bfs);

				FOR_NBORS_OF_ATOM(nbor, curr_bfs) {
					if (nbor->GetValence() == 2 && nbor->GetAtomicNum() != TREE_EXT_CONN && a->GetAtomicNum() != TREE_BRANCH_POINT && !visited_branches.HasAtom(&*nbor)) {
						to_visit.push(&*nbor);
					}
				}
			}

			branches_to_combine.push_back(a_nbors);
		}
	}
	for (std::vector<VirtualMol>::iterator branch_set=branches_to_combine.begin(); branch_set!=branches_to_combine.end(); ++branch_set) {
		// Make a new PA at the centroid of the branch sites
		OBMol branch_mol = branch_set->ToOBMol();
		vector3 branch_centroid = getCentroid(&branch_mol, false);
		PseudoAtom new_branch = formAtom(frag_molp, branch_centroid, TREE_INT_BRANCH);
		VirtualMol new_origin_map = VirtualMol(frag_map->origin_molp);

		// Keep bonds, remove original branch PA's, and update MappedMol accounting
		VirtualMol visited_nbors(frag_molp);
		AtomSet set_atoms = branch_set->GetAtoms();
		for (AtomSet::iterator it=set_atoms.begin(); it!=set_atoms.end(); ++it) {
			FOR_NBORS_OF_ATOM(n, *it) {
				if (!branch_set->HasAtom(&*n) && !visited_nbors.HasAtom(&*n)) {
					visited_nbors.AddAtom(&*n);
					formBond(frag_molp, new_branch, &*n, 1);
				}
			}
			frag_molp->DeleteAtom(*it);
			new_origin_map.AddVirtualMol(frag_map->copy_pa_to_multiple[*it]);
			frag_map->copy_pa_to_multiple.erase(*it);
		}
		frag_map->copy_pa_to_multiple[new_branch] = new_origin_map;
	}
	if (DEBUG_WITH_CIFS) {writeCIF(frag_molp, GetOutputPath("debug_tree_8_internal_branches.cif")); }


	// Assign everything else as a branch point.
	// Connections with 2+ explicit valence (and an implicit connection bond) have already been assigned as BPs.
	FOR_ATOMS_OF_MOL(a, frag_molp) {
		if (a->GetAtomicNum() == TREE_BRANCH_POINT) {
			// Already detected
		} else if (a->GetAtomicNum() != TREE_INT_BRANCH && a->GetAtomicNum() != TREE_EXT_CONN) {
			a->SetAtomicNum(TREE_BRANCH_POINT);
			if (a->GetValence() < 3) {
				obErrorLog.ThrowError(__FUNCTION__, "Found branch point with fewer than three connections.", obError);
			}
		}
	}
	if (DEBUG_WITH_CIFS) {writeCIF(frag_molp, GetOutputPath("debug_tree_9_branch_points.cif")); }
}  // end of AllNodeDeconstructor::TreeDecomposition


void AllNodeDeconstructor::CollapseRings(MappedMol *fragment_to_simplify, bool fuse_bridged_rings) {
	// Collapse rings in the MappedMol to single PA's, optionally fusing bridged rings
	// together instead of keeping them as two distinct entities

	// Convenient aliases
	OBMol* frag_molp = &(fragment_to_simplify->mol_copy);
	OBUnitCell* uc = getPeriodicLattice(frag_molp);

	// Get the list of rings to collapse
	std::vector<VirtualMol> rings_to_simplify;
	// FOR_RINGS_OF_MOL was causing weird segfaults in OBMol's without rings.
	// Unclear if it's a bug in FOR_RINGS_OF_MOL or not explicitly running ring perception.
	// At any rate, using the results of GetSSSR, instead, avoids these segfaults
	std::vector<OBRing*> ring_vec = frag_molp->GetSSSR();
	//FOR_RINGS_OF_MOL(r, *frag_molp) {
	for (std::vector<OBRing*>::iterator r=ring_vec.begin(); r!=ring_vec.end(); ++r) {
		// Get ring atoms
		VirtualMol current_ring(frag_molp);
		std::vector<int> r_atoms = (*r)->_path;
		for (std::vector<int>::iterator r_a=r_atoms.begin(); r_a!=r_atoms.end(); ++r_a) {
			current_ring.AddAtom(frag_molp->GetAtom(*r_a));
		}
		rings_to_simplify.push_back(current_ring);
	}

	// If the fuse_bridged_rings option is enabled, join together rings sharing at least two bridged atoms
	// There is probably an algorithmically faster way to accomplish the same goal, e.g. something related
	// to [heap queues](https://github.com/python/cpython/blob/3.6/Lib/heapq.py), but let's avoid making
	// premature optimizations to the code.
	if (fuse_bridged_rings) {
		int fuse_simplifications = 0;
		do {
			fuse_simplifications = 0;
			// Use integers instead of vector iterators to make the code easier to follow
			int num_rings = rings_to_simplify.size();
			for (int i = 0; i < num_rings; ++i) {
				for (int j = i+1; j < num_rings; ++j) {
					VirtualMol i_ring = rings_to_simplify[i];  // refresh in the inner loop in case it's been modified by a previous j iteration
					VirtualMol j_ring = rings_to_simplify[j];
					AtomSet j_set = j_ring.GetAtoms();
					int ij_overlap = 0;
					for (AtomSet::iterator j_it=j_set.begin(); j_it!=j_set.end(); ++j_it) {
						if (i_ring.HasAtom(*j_it)) {
							++ij_overlap;
						}
					}
					if (ij_overlap >= 2) {
						i_ring.AddVirtualMol(j_ring);
						rings_to_simplify[i] = i_ring;
						rings_to_simplify[j] = VirtualMol(frag_molp);  // erase the ring, which will avoid double counting
					}
				}
			}
			for (int k = num_rings-1; k >= 0; --k) {  // Traverse in reverse order to avoid invalidating the index
				if (rings_to_simplify[k].NumAtoms() == 0) {
					rings_to_simplify.erase(rings_to_simplify.begin() + k);
					++fuse_simplifications;
				}
			}
		} while (fuse_simplifications);
	}

	// Collapse the rings
	for (std::vector<VirtualMol>::iterator r=rings_to_simplify.begin(); r!=rings_to_simplify.end(); ++r) {
		AtomSet r_nbors;
		AtomSet r_atoms = r->GetAtoms();

		OBMol r_mol = r->ToOBMol();
		vector3 r_loc = getCentroid(&r_mol, false);

		VirtualMol r_origin_atoms(fragment_to_simplify->origin_molp);

		for (AtomSet::iterator it=r_atoms.begin(); it!=r_atoms.end(); ++it) {
			FOR_NBORS_OF_ATOM(nbor, **it) {
				if (!r->HasAtom(&*nbor)) {  // If it's not part of the ring
					r_nbors.insert(&*nbor);
				}
			}

			r_origin_atoms.AddVirtualMol(fragment_to_simplify->copy_pa_to_multiple[*it]);
			fragment_to_simplify->copy_pa_to_multiple.erase(*it);

			frag_molp->DeleteAtom(*it);
		}

		// Replace the ring atoms with PA connected to the old neighbors
		PseudoAtom new_ring_pa = formAtom(frag_molp, r_loc, TREE_PA_ELEMENT);
		fragment_to_simplify->copy_pa_to_multiple[new_ring_pa] = r_origin_atoms;
		for (AtomSet::iterator nbor_to_bond=r_nbors.begin(); nbor_to_bond!=r_nbors.end(); ++nbor_to_bond) {
			formBond(frag_molp, new_ring_pa, *nbor_to_bond, 1);
		}
	}
}


void AllNodeDeconstructor::WriteCIFs() {
	// Call base class exporter, plus new outputs for branch info
	SingleNodeDeconstructor::WriteCIFs();
	branch_points_orig.ToCIF(GetOutputPath("branch_points.cif"));
	branches_orig.ToCIF(GetOutputPath("branches.cif"));
	writeCIF(&branch_points_pa, GetOutputPath("branch_points_simplified.cif"));
	writeCIF(&branches_pa, GetOutputPath("branches_simplified.cif"));
}


} // end namespace OpenBabel
