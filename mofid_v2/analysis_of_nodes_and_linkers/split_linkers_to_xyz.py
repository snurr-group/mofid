from ase.io import read as ase_read
from ase.io import write as ase_write
from ase import neighborlist
from ase.formula import Formula
import networkx as nx
import glob
from ase.build import sort

def dict2str(dct):
    """Convert symbol-to-number dict to str.
    """
    return ''.join(symb + (str(n)) for symb, n in dct.items())

metals  = ['Li', 'Na', 'K', 'Rb', 'Cs', 'Fr', 'Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra',
          'Al', 'Ga', 'Ge', 'In', 'Sn', 'Sb', 'Tl', 'Pb', 'Bi', 'Po',
          'Sc', 'Ti', 'V' , 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
          'Y' , 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
          'Hf', 'Ta', 'W' , 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
          'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'U', 'Tm', 'Yb', 'Lu',
          'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr'
         ]

def split_nodes_from_cif(cifpath):
    try:
        atoms = ase_read(cifpath)
    except:
        print('Error with reading CIF in {}'.format(cifpath))
        return 1

    cutOff = neighborlist.natural_cutoffs(atoms)
    neighborList = neighborlist.NeighborList(cutOff, self_interaction=False, bothways=True, skin=0.3)
    neighborList.update(atoms)
    G = nx.Graph()
    
    for k in range(len(atoms)):
        tup = (k, {"element":"{}".format(atoms.get_chemical_symbols()[k]), "pos": atoms.get_positions()[k]})
        G.add_nodes_from([tup])

    for k in range(len(atoms)):
        for i in neighborList.get_neighbors(k)[0]:
            G.add_edge(k, i)
    
    Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
    
    form_dicts = []
    for index, g in enumerate(Gcc):
        g = list(g)
        fragment = atoms[g]
        fragment = sort(fragment)
        form_dict = fragment.symbols.formula.count()
        form_dicts.append(dict2str(form_dict))
    
    nodes = []
    unique_formdicts = []

    if len(form_dicts) > 1:
        for index, form_dict in enumerate(form_dicts):
            if form_dict not in unique_formdicts:
                nodes.append(atoms[list(Gcc[index])])
                unique_formdicts.append(form_dict)
    elif len(form_dicts) == 1:
        nodes.append(atoms[list(Gcc[0])])
        unique_formdicts.append(form_dicts[0])

    for index, atom in enumerate(nodes):
        print(cifpath, cifpath.split('/')[0])
        cifname = cifpath.split('/')[0]
        ase_write('{}/{}_MetalOxolinker{}.xyz'.format(cifname, cifname, index), nodes[index])
    return 0

paths = glob.glob('*/MetalOxo/linkers.cif')
for path in paths:
    print(path)
    cifname = path.split('/')[0]
    check = split_nodes_from_cif(path)
    if check == 1:
        with open('failed_refs.txt', 'a') as f:
            f.write('{}\n'.format(cifname))
