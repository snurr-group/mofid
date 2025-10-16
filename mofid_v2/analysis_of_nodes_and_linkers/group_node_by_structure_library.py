import copy
import ase
import json
from io import StringIO
from ase.io import read, write
import numpy as np
import collections
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.analysis.structure_matcher import ElementComparator
from pymatgen.analysis.structure_matcher import StructureMatcher

def convert_ase_pymat(ase_objects):
    structure_lattice = Lattice(ase_objects.cell)
    structure_species = ase_objects.get_chemical_symbols()
    structure_positions = ase_objects.get_positions()
    return Structure(structure_lattice,structure_species,structure_positions)

def remove_pbc_cuts(atoms):
    """Remove building block cuts due to periodic boundary conditions. After the
    removal, the atoms object is centered at the center of the unit cell.

    Parameters
    ----------
    atoms : ase.Atoms
        The atoms object to be processed.

    Returns
    -------
    ase.Atoms
        The processed atoms object.
    """
    try:
        # setting cuttoff parameter
        scale = 1.4
        cutoffs = ase.neighborlist.natural_cutoffs(atoms)
        cutoffs = [scale * c for c in cutoffs]

        #making neighbor_list
        #graphs of single metal is not constructed because there is no neighbor.
        I, J, D = ase.neighborlist.neighbor_list("ijD",atoms,cutoff=cutoffs)

        nl = [[] for _ in atoms]
        for i, j, d in zip(I, J, D):
            nl[i].append((j, d))
        visited = [False for _ in atoms]
        q = collections.deque()

        # Center of the unit cell.
        abc_half = np.sum(atoms.get_cell(), axis=0) * 0.5

        positions = {}
        q.append((0, np.array([0.0, 0.0, 0.0])))
        while q:
            i, pos = q.pop()
            visited[i] = True
            positions[i] = pos
            for j, d in nl[i]:
                if not visited[j]:
                    q.append((j, pos + d))
                    visited[j] = True

        centroid = np.array([0.0, 0.0, 0.0])
        for v in positions.values():
            centroid += v
        centroid /= len(positions)

        syms = [None for _ in atoms]
        poss = [None for _ in atoms]

        for i in range(len(atoms)):
            syms[i] = atoms.symbols[i]
            poss[i] = positions[i] - centroid + abc_half


        atoms = ase.Atoms(
            symbols=syms, positions=poss, pbc=True, cell=atoms.get_cell()
        )
        
        # resize of cell
        cell_x = np.max(atoms.positions[:,0]) - np.min(atoms.positions[:,0])
        cell_y = np.max(atoms.positions[:,1]) - np.min(atoms.positions[:,1])
        cell_z = np.max(atoms.positions[:,2]) - np.min(atoms.positions[:,2])
        cell = max([cell_x,cell_y,cell_z])

        atoms.set_cell([cell+2,cell+2,cell+2, 90,90,90])
        center_mass = atoms.get_center_of_mass()
        cell_half  = atoms.cell.cellpar()[0:3]/2
        atoms.positions = atoms.positions - center_mass + cell_half    

        return atoms
    except:
        return atoms

def return_xyz_list(atoms):
    xyz_stringio = StringIO()
    write(xyz_stringio, atoms, format='xyz')
    xyz_stringio.seek(0)
    return xyz_stringio.readlines()


if __name__ == "__main__":
    f = open("nodes_details.json")
    node_jsons = json.load(f)


    type_list_in_formula = {}
    for formula in node_jsons.keys():

        node_list = sorted(node_jsons[formula])
        unique_list_in_formula = []
        if len(node_list) >1:

            i = 1
            type_list = {}
            while len(node_list) > 0:

                duplicate_list_in_formula = [node_list[0]]
                node_list_for_loop = copy.copy(node_list)
                mof = remove_pbc_cuts(read(node_list[0] + ".xyz"))
                write(formula+"_Type-" + str(i)+".xyz",mof)
                a = convert_ase_pymat(mof)
                node_list_for_loop.pop(0)

                for node in reversed(node_list_for_loop):
                    b = convert_ase_pymat(remove_pbc_cuts(read(node + ".xyz")))
                    matcher = StructureMatcher(ltol = 0.3, stol = 2, angle_tol = 5, primitive_cell = False, scale = False, comparator=ElementComparator())
                    are_structures_matching = matcher.fit(a, b)
                    if are_structures_matching:
                        node_list_for_loop.pop(node_list_for_loop.index(node))
                        duplicate_list_in_formula.append(node)

                type_list["Type-" + str(i)] = {"node-file":sorted(duplicate_list_in_formula), "xyz": return_xyz_list(mof)}
                node_list = node_list_for_loop
                i += 1
            type_list_in_formula[formula] = type_list

        elif len(node_list) ==1:
            mof = remove_pbc_cuts(read(node_list[0] + ".xyz"))
            write(formula+"_Type-1.xyz",mof)
            type_list_in_formula[formula] = {"Type-1":{"node-file":node_list[0], "xyz": return_xyz_list(mof)}} #

    with open("node_grouped_by_structure.json", 'w') as file:
        json.dump(type_list_in_formula, file,indent=2)