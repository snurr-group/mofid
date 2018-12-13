#include "framework.h"
#include "obdetails.h"
#include "periodic.h"

#include <string>
#include <vector>
#include <queue>

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/generic.h>
#include <openbabel/obiter.h>
#include <openbabel/elements.h>
#include <openbabel/obconversion.h>
#include <openbabel/phmodel.h>


namespace OpenBabel
{

bool importCIF(OBMol* molp, std::string filepath, bool bond_orders, bool makeP1) {
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

	// Strip all of the original CIF labels, so they don't interfere with the automatically generated labels in the output
	FOR_ATOMS_OF_MOL(a, *molp) {
		if (a->HasData("_atom_site_label")) {
			a->DeleteData("_atom_site_label");
		}
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

OBMol initMOFwithUC(OBMol *orig_in_uc) {
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
		OBMol candidate = initMOFwithUC(mol);  // Have to check the match for infinite rods
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

} // end namespace OpenBabel
