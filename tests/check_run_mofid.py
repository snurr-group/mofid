import os
from mofid.run_mofid import cif2mofid
from mofid.paths import resources_path

test_mofs = os.path.join(resources_path,'TestCIFs')
mof1 = os.listdir(test_mofs)[0]
cif_path = os.path.join(test_mofs,mof1)
mofid = cif2mofid(cif_path)