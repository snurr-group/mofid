import os
from mofid.run_mofid import cif2mofid
from mofid.paths import resources_path

cif_path = os.path.join(resources_path,'TestCIFs','ABAVIJ_clean.cif')
mofid = cif2mofid(cif_path)
if (not mofid['mofkey'].startswith('Co.TWBYWOBDOCUKOW.MOFkey-v1.rtl')
    or mofid['mofid'].startswith('[Co].[O-]C(=O)c1ccncc1 MOFid-v1.rtl.cat0;ABAVIJ_clean')):
    raise ValueError('Test failed!')
