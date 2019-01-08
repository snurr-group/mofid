// Main MOFid code to decompose a CIF into nodes, linkers, and solvents.
// Writes the SMILES and catenation info to stdout, errors to stderr, and
// relevant CIFs and simplified topology.cgd to the Test/ directory.

// See https://openbabel.org/docs/dev/UseTheLibrary/CppExamples.html
// Get iterator help from http://openbabel.org/dev-api/group__main.shtml
// Visualize with http://baoilleach.webfactional.com/site_media/blog/emscripten/openbabel/webdepict.html

// Instead of including code to display the full MOF SMILES, just use openbabel natively:
// obabel CIFFILE -ap -ocan

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <string>

#include <openbabel/mol.h>
#include <openbabel/obiter.h>
#include <openbabel/babelconfig.h>
#include <openbabel/elements.h>

#include "config_sbu.h"
#include "invector.h"
#include "obdetails.h"
#include "deconstructor.h"
#include "framework.h"
#include "periodic.h"
#include "pseudo_atom.h"
#include "topology.h"
#include "virtual_mol.h"


using namespace OpenBabel;  // See http://openbabel.org/dev-api/namespaceOpenBabel.shtml


// Function prototypes
std::string analyzeMOF(std::string filename);
extern "C" void analyzeMOFc(const char *cifdata, char *analysis, int buflen);
extern "C" int SmilesToSVG(const char* smiles, int options, void* mbuf, unsigned int buflen);


int main(int argc, char* argv[])
{
	obErrorLog.SetOutputLevel(obWarning);  // See also http://openbabel.org/wiki/Errors
	char* filename;
	filename = argv[1];  // TODO: Check usage later

	// Set up the babel data directory to use a local copy customized for MOFs
	// (instead of system-wide Open Babel data)
	std::stringstream dataMsg;
	dataMsg << "Using local Open Babel data saved in " << LOCAL_OB_DATADIR << std::endl;
	obErrorLog.ThrowError(__FUNCTION__, dataMsg.str(), obAuditMsg);
	// Use setenv instead of putenv, per advice about string copies vs. pointers: http://stackoverflow.com/questions/5873029/questions-about-putenv-and-setenv/5876818#5876818
	// This is similar to the approach of cryos/avogadro:main.cpp:127
	// Per my objective, this only sets the environment within the scope of the sbu.exe program.
	// But Windows defines a separate _putenv, etc: https://stackoverflow.com/questions/17258029/c-setenv-undefined-identifier-in-visual-studio
#ifdef _WIN32
	_putenv_s("BABEL_DATADIR", LOCAL_OB_DATADIR);
#else
#ifndef __INTELLISENSE__ // Ignore setenv error in vscode:
	setenv("BABEL_DATADIR", LOCAL_OB_DATADIR, 1);
#endif  // vscode error workaround
#endif

	std::string mof_results = analyzeMOF(std::string(filename));
	if (mof_results == "") {  // No MOFs found
		return(1);
	} else {
		std::cout << mof_results;
		return(0);
	}
}

std::string analyzeMOF(std::string filename) {
	// Extract components of the MOFid
	// Reports nodes/linkers, number of nets found, and writes CIFs to Test/ directory

	OBMol orig_mol;
	// Massively improving performance by skipping kekulization of the full MOF
	if (!importCIF(&orig_mol, filename, false)) {
		std::cerr << "Error reading file: %s" << filename << std::endl;
		return "";
	}

	// TODO: Write out the original mol like this (near instantly) for debugging (maybe as part of importCIF)
	writeCIF(&orig_mol, "Test/orig_mol.cif");
	std::ofstream file_info;
	file_info.open("Test/mol_name.txt", std::ios::out | std::ios::trunc);
	if (file_info.is_open()) {
		file_info << filename << std::endl;
		file_info.close();
	}

	MOFidDeconstructor simplifier(&orig_mol);
	simplifier.SetOutputDir("Test/");
	simplifier.SimplifyMOF();
	simplifier.WriteCIFs();

	return simplifier.GetMOFInfo();
}

extern "C" {
void analyzeMOFc(const char *cifdata, char *analysis, int buflen) {
	// Wrap analyzeMOF with C compatibility for Emscripten usage
	// Make the return value an input param, since c_str is (likely?) a pointer to
	// an internal data structure in std::string, which will go out of scope.
	// A good buflen for the output might be 2^16, or 65536
	const char *TEMP_FILE = "from_emscripten.cif";
	std::ofstream cifp(TEMP_FILE, std::ios::out | std::ios::trunc);
	cifp << std::string(cifdata);
	cifp.close();

	strncpy(analysis, analyzeMOF(TEMP_FILE).c_str(), buflen);
}

int SmilesToSVG(const char* smiles, int options, void* mbuf, unsigned int buflen) {
	// From Noel O'Boyle: https://baoilleach.blogspot.com/2015/02/cheminformaticsjs-open-babel.html
	// See also http://baoilleach.webfactional.com/site_media/blog/emscripten/openbabel/webdepict.html
	OpenBabel::OBMol mol;
	OpenBabel::OBConversion conv;
	conv.SetInFormat("smi");

	bool ok = conv.ReadString(&mol, smiles);
	if (!ok) return -1;

	if (options==0) {
		conv.SetOutFormat("ascii");
		conv.AddOption("a", OpenBabel::OBConversion::OUTOPTIONS, "2.2");

	} else {
		conv.SetOutFormat("svg");
		conv.AddOption("C", OpenBabel::OBConversion::OUTOPTIONS);
		conv.AddOption("P", OpenBabel::OBConversion::OUTOPTIONS, "500");
	}

	std::string out = conv.WriteString(&mol);

	if (out.size()+1 >= buflen)
			return -1;

	char* dst = (char*)mbuf;
	strcpy(dst, out.c_str());
	return out.size();
}
}  // extern "C"

