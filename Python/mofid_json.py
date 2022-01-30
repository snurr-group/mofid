"""
Adapting run_mofid.py to return the full json output

Returning the cif2mofid output from within a docker container back to the host
without requiring a full compiler environment and Python package install.
See README for usage details.

@author: Ben Bucior and Andrew Rosen
"""

import sys
import json
import tempfile
from mofid.run_mofid import cif2mofid, DEFAULT_OUTPUT_PATH


if __name__ == '__main__':
	output_path = DEFAULT_OUTPUT_PATH
	args = sys.argv[1:]
	with tempfile.NamedTemporaryFile('w') as fp:
		if len(args) == 0:  # reading from stdin
			for line in sys.stdin.read():
				fp.write(line)
			fp.seek(0)
			cif_file = fp.name
			# This code now works, though the filename tag (thus MOFid) is messed up
			# Could consider an assert for the ;filename and then replace it with something
		else:  # same format as before: python mofid_json.py path_to_cif_for_analysis.cif OutputPathIfNonstandard
			cif_file = args[0]
			if len(args) == 2:
				output_path = args[1]

		identifiers = cif2mofid(cif_file, output_path)
	print(json.dumps(identifiers))
