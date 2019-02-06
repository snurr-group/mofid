#include "periodic.h"

#include <map>
#include <queue>

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/obiter.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/generic.h>


namespace OpenBabel
{

OBUnitCell* getPeriodicLattice(OBMol *mol) {
	// Replacement for the old OBMol.GetPeriodicLattice helper function
	return (OBUnitCell*)mol->GetData(OBGenericDataType::UnitCell);
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

int3 GetPeriodicDirection(OBBond *bond) {
	// What is the unit cell of the end atom wrt the first?
	// Returns {0,0,0} if not periodic or if wrapping is not required.

	int3 direction(0, 0, 0);

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

UCMap unwrapFragmentUC(OBMol *fragment, bool allow_rod, bool warn_rod) {
	// Starting with a random atom, traverse the fragment atom-by-atom to determine
	// which unit cell each atom belongs to in a self-consistent manner.
	// Includes an optional parameter to allow 1D periodic fragments (e.g. MIL-47).
	// By default, these are forbidden (since the ordering is undefined) and returns an empty map.

	std::queue<OBAtom*> to_visit;
	UCMap unit_cells;
	if (fragment->NumAtoms() == 0) {
		return unit_cells;
	}

	// Start at whichever atom is (randomly?) saved first
	// Note: atom arrays begin with 1 in OpenBabel, while bond arrays begin with 0.
	OBAtom* start_atom = fragment->GetAtom(1);
	to_visit.push(start_atom);
	unit_cells[start_atom] = int3(0, 0, 0);  // original unit cell

	while (!to_visit.empty()) {
		OBAtom* current = to_visit.front();
		to_visit.pop();
		FOR_NBORS_OF_ATOM(nbr, current) {
			OBBond* nbr_bond = fragment->GetBond(current, &*nbr);
			int3 uc = GetPeriodicDirection(nbr_bond);  // TODO: Is there a good way to remove this function?  Otherwise, add it as a MOF helper method
			if (nbr_bond->GetBeginAtom() == &*nbr) {  // opposite bond direction as expected
				uc = int3(-1*uc.x, -1*uc.y, -1*uc.z);
			}

			int3 current_uc = unit_cells[current];
			uc = int3(current_uc.x + uc.x, current_uc.y + uc.y, current_uc.z + uc.z);

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
		int3 uc_shift = it->second;

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
		int3 uc_shift = it->second;  // <first: second> = <key: value> of a map/dict.
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

} // end namespace OpenBabel
