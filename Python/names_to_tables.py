"""
Break a MOFid file into tables of its components

Converts a database of MOFid strings into a tables of the detected components
where the primary key in all the tables is the name.

@author: Ben Bucior
"""

import os, sys
import copy
from extract_moffles import parse_moffles


def dict_to_delim(to_export, filename, delim="\t"):
	with open(filename, "w") as f:
		for key in to_export:
			if type(to_export[key]) is list:
				for j in to_export[key]:
					f.write(key + delim + str(j) + "\n")
			else:
				f.write(key + delim + str(to_export[key]) + "\n")


class MOFExporter:
	# Exports a .smi-formatted list of MOFid into separate tables for the various parts
	def __init__(self):
		self.tables = dict()
		self.datatypes = list()  # smiles, topology, catenation, etc.

	def parse(self, filename):
		with open(filename, 'r') as f:
			for line in f:
				parsed = parse_moffles(line)
				name = parsed['name']
				del parsed['name']
				parsed['smiles_part'] = parsed['smiles'].split('.')
				parsed['base_topology'] = parsed['topology'].split(",")[0]
				parsed['extra_topology'] = ",".join(parsed['topology'].split(",")[1:])  # '' if empty
				self.tables[name] = copy.deepcopy(parsed)
				self.datatypes = parsed.keys()
		return self

	def _tidy_tables(self):
		tidy_dict = dict()
		for info in self.datatypes:
			tidy_dict[info] = dict()

		for mof in self.tables:
			for info in self.datatypes:
				tidy_dict[info][mof] = self.tables[mof][info]

		return tidy_dict

	def write(self, folder='.'):
		if not os.path.isdir(folder):
			os.mkdir(folder)

		tidy_output = self._tidy_tables()
		for key in tidy_output:
			dict_to_delim(tidy_output[key], folder + "/" + key + ".tsv", delim="\t")



if __name__ == "__main__":
	args = sys.argv[1:]
	if len(args) != 1:
		raise SyntaxError("Extract info from a list of MOFid strings.  Only a single .smi filename expected.")

	input_smi = args[0]
	MOFExporter().parse(input_smi).write('TableOutput')
