/* tsfm_smiles: run a SMARTS reaction transformation on a SMILES string using OBChemTsfm */
/* Derived from cheminformatics.py:openbabel_replace.  See also current cpp_cheminformatics.py */

#include <sstream>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/kekulize.h>
#include <openbabel/babelconfig.h>
#include <openbabel/phmodel.h>
#include "config_sbu.h"


using namespace OpenBabel;  // See http://openbabel.org/dev-api/namespaceOpenBabel.shtml


int main(int argc, char* argv[])
{
	obErrorLog.SetOutputLevel(obInfo);  // See also http://openbabel.org/wiki/Errors
	std::string mol_smiles(argv[1]);
	std::string query(argv[2]);
	std::string replacement(argv[3]);

	// Set up the babel data directory
	std::stringstream dataMsg;
	dataMsg << "Using local Open Babel data saved in " << LOCAL_OB_DATADIR << std::endl;
	obErrorLog.ThrowError(__FUNCTION__, dataMsg.str(), obAuditMsg);
#ifdef _WIN32
	_putenv_s("BABEL_DATADIR", LOCAL_OB_DATADIR);
#else
	setenv("BABEL_DATADIR", LOCAL_OB_DATADIR, 1);
#endif

	// Read input SMILES
	OBMol orig_mol;
	OBConversion reader;
	reader.SetInFormat("smi");
	reader.AddOption("s", OBConversion::INOPTIONS);
	if (!reader.ReadString(&orig_mol, mol_smiles)) {
		obErrorLog.ThrowError(__FUNCTION__, "Error reading input SMILES", obError);
		exit(1);
	}

	// Run the transformation
	OBChemTsfm tsfm;
	if (!tsfm.Init(query, replacement)) {
		obErrorLog.ThrowError(__FUNCTION__, "Internal error: could not parse reaction transform", obError);
		exit(2);
	}
	tsfm.Apply(orig_mol);
	orig_mol.UnsetAromaticPerceived();
	OBKekulize(&orig_mol);

	// Write the output SMILES
	OBConversion writer;
	writer.SetOutFormat("can");
	writer.AddOption("i");  // Ignore SMILES chirality for now
	std::cout << writer.WriteString(&orig_mol) << std::endl;
}
