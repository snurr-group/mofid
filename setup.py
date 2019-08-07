from setuptools import setup

setup(name='mofid',
      description='A system for rapid identification and analysis of metal-organic frameworks',
      author='Benjamin J. Bucior',
      url='https://github.com/snurr-group/mofid', 
      version='1.0.0',
      packages=['mofid',],
      package_dir = {'mofid':'Python'},
      license='GNU',
      install_requires=['subprocess32>="3.5.0";python_version<"3.0"']
     )
