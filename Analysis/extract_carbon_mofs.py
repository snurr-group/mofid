# Only include data for carbon-containing MOFs

import os
import sys
import shutil
import glob
import re

DEST_DIR = 'Summary/OrigCore/'

if not os.path.isdir(DEST_DIR):
    print('Moving original CoRE MOF files to a new directory: ' + DEST_DIR)
    os.mkdir(DEST_DIR)
    core_files = glob.glob('Summary/core*')
    for core_output in core_files:
        shutil.move(core_output, DEST_DIR)

def read_file(data_file):
    with open(os.path.join('Summary/OrigCore/', data_file), 'r') as f:
        data_contents = f.readlines()
        # Let's keep the newlines so we don't have to rewrite them later:
        #data_contents = [x.strip('\n') for x in data_contents]
    return data_contents

mofids = read_file('core_mofid.smi')
#smis = read_file('core.smi')
mofkeys = read_file('core_mofkey.tsv')
elements = read_file('core_molec_formula.tsv')

out_mofids = open('Summary/core_mofid.smi', 'w')
out_mofkeys = open('Summary/core_mofkey.tsv', 'w')
out_files = open('Summary/core_extracted_files_with_carbon.txt', 'w')

# Remove the .tsv headers and redirect as necessary
elements.pop(0)
out_mofkeys.write(mofkeys.pop(0))

# Map the input to output files line-by-line
for lineno in range(len(elements)):
    cif_name = elements[lineno].split()[0]
    assert cif_name.endswith('.cif') or cif_name.endswith('.CIF')
    cif_name = cif_name[:-4]
    # Check that our lines are properly aligned between all of the MOFid code output files
    assert mofids[lineno].endswith(';' + cif_name + '\n')
    assert mofkeys[lineno].startswith(cif_name + '.cif\t')

    contains_carbon = re.search(r'[ \t]C ', elements[lineno])
    if contains_carbon:
        out_mofids.write(mofids[lineno])
        out_mofkeys.write(mofkeys[lineno])
        out_files.write(cif_name + '\n')

out_mofids.close()
out_mofkeys.close()
out_files.close()
shutil.copy2('Summary/core_mofid.smi', 'Summary/core.smi')

