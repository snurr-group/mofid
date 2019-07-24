import os
from mofid.run_mofid import cif2mofid
from mofid.paths import resources_path

cif_path = os.path.join(resources_path,'ABAVIJ_clean.cif')
mofid = cif2mofid(cif_path)
