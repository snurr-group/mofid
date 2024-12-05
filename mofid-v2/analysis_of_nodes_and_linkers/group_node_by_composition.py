from ase.io import read as ase_read
from ase.formula import Formula
import networkx as nx
import glob
from ase.build import sort
import json

def dict2str(dct):
    """Convert symbol-to-number dict to str.
    """
    return ''.join(symb + (str(n)) for symb, n in dct.items())

dict_formulas = {}

xyzfiles = sorted(glob.glob('*.xyz'))

for count, xyzfile in enumerate(xyzfiles):
    ref = xyzfile[:-4]
    atoms = ase_read(xyzfile)
    atoms = sort(atoms)
    form_dict = atoms.symbols.formula.count()
    form = dict2str(form_dict)

    if form in list(dict_formulas):
        dict_formulas[form].append(ref)
    else:
        dict_formulas[form] = [ref]
    """ 
    if count% 1000 == 0:
        print(dict_formulas)
        print(len(dict_formulas))
    """
with open('nodes_details.json', 'w') as f:
    json.dump(dict_formulas, f, indent=4)
