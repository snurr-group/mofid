import subprocess
from openbabel import openbabel as ob
import json

def dict2str(dct):
    """Convert symbol-to-number dict to str.
    """
    return ''.join(symb + (str(n)) for symb, n in dct.items())

def are_identical_smiles(smiles1, smiles2):
    obConversion = ob.OBConversion()
    obConversion.SetInFormat("smi")
    obConversion.SetOutFormat("can")  # Set output to canonical SMILES

    mol1 = ob.OBMol()
    mol2 = ob.OBMol()

    # Read the SMILES strings into molecules
    obConversion.ReadString(mol1, smiles1)
    obConversion.ReadString(mol2, smiles2)

    # Convert molecules to canonical SMILES
    can_smiles1 = obConversion.WriteString(mol1).strip()
    can_smiles2 = obConversion.WriteString(mol2).strip()

    # Compare canonical SMILES
    return can_smiles1 == can_smiles2

metals  = ['Li', 'Na', 'K', 'Rb', 'Cs', 'Fr', 'Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra',
          'Al', 'Ga', 'Ge', 'In', 'Sn', 'Sb', 'Tl', 'Pb', 'Bi', 'Po',
          'Sc', 'Ti', 'V' , 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
          'Y' , 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
          'Hf', 'Ta', 'W' , 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
          'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'U', 'Tm', 'Yb', 'Lu',
          'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr'
         ]

linkers = []
refcodes = []
dict_data = {}

with open('all_mofid-v1.txt', 'r') as f:
    for idx, line in enumerate(f):
        if idx % 100 == 0:
            print(idx, len(dict_data), dict_data)
        if line.startswith('*'):
            continue
        mofid = line.strip()
        main = mofid.split()[0]
        supp = mofid.split()[1]
        refcode = supp.split(';')[-1]
        
        main = main.split('.')
        for new_smile in main:
            check = False
            for j in metals:
                if j in new_smile:
                    check = True
                    break
            if check == True:
                continue
            else:
                check_match = False
                if idx == 0:
                    dict_data[new_smile] = [refcode]
                    continue
                for smile in list(dict_data):
                    if are_identical_smiles(new_smile, smile):
                        dict_data[smile].append(refcode)
                        check_match = True
                        break
                if check_match ==  False:
                    dict_data[new_smile] = [refcode]

with open('smiles_by_ref.json', 'w') as f:
    json.dump(dict_data, f, indent=4)
