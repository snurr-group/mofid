from setuptools import setup

setup(name='mofid',
      description='A system for rapid identification and analysis of metal-organic frameworks',
      author='Benjamin J. Bucior',
      url='https://github.com/swanickt/mofid',
      version='1.1.1',
      packages=['mofid', 'mofid_v2', 'mofid_v2.analysis_of_nodes_and_linkers'],
      package_dir = {'mofid':'Python'},
      license='GNU',
      install_requires=[
            'numpy',
            'ase',
            'networkx',
            'pymatgen',
      ],
      python_requires='>=3.7',
      #install_requires=['subprocess32>="3.5.0";python_version<"3.0"']
     )
