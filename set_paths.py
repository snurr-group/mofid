import os

mofid_path = os.path.split(os.path.abspath(__file__))[0]
mofidpy_path = os.path.join(mofid_path,'Python')
bin_path = os.path.join(mofid_path,'bin')
openbabel_path = os.path.join(mofid_path,'openbabel')
resources_path = os.path.join(mofid_path,'Resources')

for p in [mofid_path,mofidpy_path,bin_path,openbabel_path,resources_path]:
    if not os.path.exists(p):
        raise ValueError('No directory '+p)

with open(os.path.join(mofidpy_path,'paths.py'),'w') as f:
    f.write("mofid_path = r'"+mofid_path+"'\nbin_path = r'"+bin_path+"'\nopenbabel_path = r'"+openbabel_path+"'\nresources_path = r'"+resources_path+"'\n")
