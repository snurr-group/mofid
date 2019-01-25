#include "deconstructor.h"
#include "invector.h"
#include "obdetails.h"
#include "framework.h"
#include "periodic.h"
#include "topology.h"

#include <string>
#include <sstream>
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

	if (write_intermediate_cifs) { WriteSimplifiedNet("simplified_test_orig.cif"); }
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
}



SingleNodeDeconstructor::SingleNodeDeconstructor(OBMol* orig_mof) : Deconstructor(orig_mof) {
	points_of_extension = VirtualMol(orig_mof);
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

	// The second VirtualMol has point(s) of extension.  The first VirtualMol is the ring minus those points.
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
				ring.RemoveAtom(curr_ring_atom);
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
			// Check coordination environment to classify as part of the metal oxide or "other"
			VirtualMol attached_hydrogens(nn->GetParent());
			OBAtom* attached_nonmetal = NULL;
			FOR_NBORS_OF_ATOM(n, *nn) {
				if (!nodes.HasAtom(&*n)) {  // metals and other bound oxygens
					if (n->GetAtomicNum() == 1) {
						attached_hydrogens.AddAtom(&*n);
					} else if (n->GetAtomicNum() == 8) {
						obErrorLog.ThrowError(__FUNCTION__, "Found a loosely bound peroxide, which will be handled by 1-c solvent codes.", obWarning);
					} else {
						if (attached_nonmetal) {  // O is already bound to a different nonmetal!
							obErrorLog.ThrowError(__FUNCTION__, "Found a metal-bound oxygen with multiple nonmetal neighbors!  Results may have inconsistencies.", obError);
						}
						attached_nonmetal = &*n;
					}
				}
			}

			if (!attached_nonmetal) {
				// Keep the oxygen and its bound H's
				nodes.AddVirtualMol(attached_hydrogens);
			} else {
				// Test the oxygen binding mode to classify it as part of node vs. linker
				OBAtom* oxygen_bridge = NULL;
				FOR_NBORS_OF_ATOM(n, *attached_nonmetal) {
					if ((temp_node.find(&*n) != temp_node.end()) && (nn != &*n)) {
						// Use temp_node because it's the original list of metals and NN's, not added carbons, rings, etc.
						// Avoids a bug when handling self-connected paddlewheels, e.g. hMOF-22242.
						// Note: the carboxylates will become "node bridges" in this instance to avoid a self-connected, periodic node.
						if (n->GetAtomicNum() != 8) {
							obErrorLog.ThrowError(__FUNCTION__, "Unexpectedly found an M-O-NM-NM-M bridge.  Excepted a second oxygen.", obError);
						} else {
							if (oxygen_bridge) {  // already defined
								obErrorLog.ThrowError(__FUNCTION__, "Unexpectedly found multiple O2's for M-O1-NM-O2-M.", obWarning);
							}
							oxygen_bridge = &*n;
						}
					}
				}
				if (!oxygen_bridge) {
					// Binding end-on, e.g. carboxylates with a single oxygen binding or catechols.
					// TODO: think about cases where an oxygen binds to multiple metal atoms.
					// Should the oxygen in that case be part of the node or linker?
					nodes.RemoveAtom(nn);  // not part of an SBU cluster
				} else {
					// If there is a bridge, include the carboxylate carbon as part of the SBU
					points_of_extension.AddAtom(attached_nonmetal);
					visited_bridges.AddAtom(nn);
					visited_bridges.AddAtom(oxygen_bridge);
					// Get hydrogens attached to the second oxygen atom
					FOR_NBORS_OF_ATOM(b_h, *oxygen_bridge) {
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
					if (bridge_bw_metals) {  // already defined
						obErrorLog.ThrowError(__FUNCTION__, "Unexpectedly found N-N-N species in node.", obWarning);
					}
					bridge_bw_metals = &*n;
				}
			}
			if (!bridge_bw_metals) {  // Binding end-on, e.g. pillared MOFs, porphyrins, and probably amines
				nodes.RemoveAtom(nn);  // not part of an SBU cluster
			} else {  // If there is a bridge, grab the whole ring
				std::pair<VirtualMol,VirtualMol> ring_info = CalculateNonmetalRing(nn, bridge_bw_metals);
				nodes.AddVirtualMol(ring_info.first);
				points_of_extension.AddVirtualMol(ring_info.second);
				visited_bridges.AddAtom(nn);
				visited_bridges.AddAtom(bridge_bw_metals);
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
	WriteSBUs("node_sbus_no_ext_conns.cif", false);
	WriteSBUs("node_sbus_with_ext_conns.cif", true);
}


void SingleNodeDeconstructor::WriteSBUs(const std::string &base_filename, bool external_bond_pa) {
	// Write out the combined SBU's, including points of extension.
	// If external_bond_pa is true, also write out pseudoatoms corresponding to the next
	// external atom connected to the point of extension.

	// Combine SBU's from nodes and PoE's (assigned as linkers for topological convenience)
	VirtualMol nodes = simplified_net.PseudoToOrig(simplified_net.GetAtomsOfRole("node"));
	VirtualMol sbus = nodes;
	sbus.AddVirtualMol(points_of_extension);
	OBMol sbu_mol = sbus.ToOBMol();  // Copy POE's as new atoms

	// Map PoE's from the original MOF to sbu_mol
	std::map<OBAtom*, OBAtom*> v_to_mol_poe;
	AtomSet poe_set = points_of_extension.GetAtoms();
	for (AtomSet::iterator it=poe_set.begin(); it!=poe_set.end(); ++it) {
		v_to_mol_poe[*it] = atomInOtherMol(*it, &sbu_mol);
	}

	// Find PoE-PoE bonds then delete them
	FOR_BONDS_OF_MOL(b, *parent_molp) {
		OBAtom* begin = b->GetBeginAtom();
		OBAtom* end = b->GetEndAtom();
		if (points_of_extension.HasAtom(begin) && points_of_extension.HasAtom(end)) {
			OBAtom* mol_begin = v_to_mol_poe[begin];
			OBAtom* mol_end = v_to_mol_poe[end];
			sbu_mol.DeleteBond(sbu_mol.GetBond(mol_begin, mol_end));
		}
	}

	if (external_bond_pa) {
		OBUnitCell* uc = getPeriodicLattice(parent_molp);
		for (AtomSet::iterator it=poe_set.begin(); it!=poe_set.end(); ++it) {
			OBAtom* poe_orig_atom = *it;  // PoE in the original MOF OBMol
			FOR_NBORS_OF_ATOM(nbor, *poe_orig_atom) {
				if (!nodes.HasAtom(&*nbor)) {  // external
					// Calculate a distance 1/3 of the way from the PoE to the external atom to avoid overlapping atoms
					// when two POEs are next to each other in the MOF (e.g. POXHER_clean.cif)
					vector3 poe_pos = poe_orig_atom->GetVector();
					vector3 external_pos = uc->UnwrapCartesianNear(nbor->GetVector(), poe_pos);
					vector3 conn_pos = uc->WrapCartesianCoordinate((2.0*poe_pos + external_pos) / 3.0);

					PseudoAtom new_conn_atom = formAtom(&sbu_mol, conn_pos, POE_EXTERNAL_ELEMENT);
					OBAtom* poe_mol_atom = v_to_mol_poe[poe_orig_atom];
					formBond(&sbu_mol, poe_mol_atom, new_conn_atom, 1);
				}
			}
		}
	}

	writeCIF(&sbu_mol, GetOutputPath(base_filename));
}



AllNodeDeconstructor::AllNodeDeconstructor(OBMol* orig_mof) : SingleNodeDeconstructor(orig_mof) {
	branches = VirtualMol(orig_mof);
	branch_points = VirtualMol(orig_mof);
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
		for (std::map<PseudoAtom,VirtualMol>::iterator it=frag_all_node.copy_pa_to_multiple.begin(); it!=frag_all_node.copy_pa_to_multiple.end(); ++it) {
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
				branch_points.AddVirtualMol(simplified_net.PseudoToOrig(collapsed));  // TODO: also consider the simplified atoms instead of original

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

				VirtualMol pa_network = frag_all_node.copy_pa_to_multiple[&*pa];
				PseudoAtom collapsed = simplified_net.CollapseFragment(pa_network);
				simplified_net.SetRoleToAtom("linker", collapsed);
				branches.AddVirtualMol(simplified_net.PseudoToOrig(collapsed));

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
		for (AtomSet::iterator it=orig_from_pa.begin(); it!=orig_from_pa.end(); ++it) {
			if (connection_points.HasAtom(*it)) {
				pa->SetAtomicNum(TREE_EXT_CONN);
			}
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
			if (a->GetValence() == 1 && a->GetAtomicNum() != TREE_EXT_CONN) {
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
	branch_points.ToCIF(GetOutputPath("branch_points.cif"));
	branches.ToCIF(GetOutputPath("branches.cif"));
}


} // end namespace OpenBabel
