/**********************************************************************
deconstructor.h - Deconstruct a MOF into its building blocks and simplified net
***********************************************************************/

#ifndef DECONSTRUCTOR_H
#define DECONSTRUCTOR_H

#include <string>
#include <utility>  // std::pair

#include <openbabel/babelconfig.h>
#include <openbabel/generic.h>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>

#include "topology.h"

namespace OpenBabel
{
// forward declarations
class OBMol;


// Constants:
// Default directory for CIF/Systre outputs.  Also used in Python/run_mofid.py and bin/sbu
const std::string DEFAULT_OUTPUT_PATH = "Output/";
const std::string NO_SBU_SUFFIX = "/NoSBU";  // base MOFidDeconstructor class
const std::string SINGLE_NODE_SUFFIX = "/SingleNode";
const std::string ALL_NODE_SUFFIX = "/AllNode";

// Default placeholder topology and details for MOFkey
const std::string DEFAULT_MOFKEY_TOPOLOGY = "";  // Alternatively, "OPTIONAL_TOPOLOGY" for user-friendliness // TODO: MAYBE NA???
const std::string MOFKEY_VERSION = "v1";
const std::string MOFKEY_SEP = ".";
const std::string MOFKEY_METAL_DELIM = ",";  // If multiple metal nodes, the delimiter between elements
const std::string MOFKEY_NO_METALS = "NA";
const std::string MOFKEY_NO_LINKERS = "MISSING_LINKERS";

// PA's to describe points of extension
const int POE_EXTERNAL_ELEMENT = 118;  // Og
const int SBU_EXTERNAL_ELEMENT = 117;  // Ts

// Atom type codes used by the MappedMol in AllNodeDeconstructor::TreeDecomposition()
const int TREE_PA_ELEMENT = 118;  // Og
const int TREE_INT_BRANCH = 117;  // Ts
const int TREE_BRANCH_POINT = 116;  // Lv
const int TREE_EXT_CONN = 115;  // Mc


// Function prototypes
std::string writeFragments(std::vector<OBMol> fragments, OBConversion obconv);
std::string exportNormalizedMol(OBMol fragment, OBConversion obconv);
std::string getSMILES(OBMol fragment, OBConversion obconv);


class Deconstructor {
// Base class for MOF deconstruction algorithms, to go from an OBMol of original atoms
// to a simplified net, its topology, and the mapping of net pseudoatoms back to the MOF.
private:  // hidden from derived classes, too
	std::string output_dir;

	Deconstructor(const Deconstructor& other);  // again, remove copy capabilities to avoid implementing them
	Deconstructor& operator=(const Deconstructor&);

protected:
	OBMol *parent_molp;
	Topology simplified_net;
	OBConversion obconv;
	bool infinite_node_detected;

	VirtualMol points_of_extension;  // Placeholder to track SBU points of extension separately from atom roles

	virtual void InitOutputFormat();
	static std::string GetBasicSMILES(OBMol fragment);

	// Network simplification steps
	virtual void DetectInitialNodesAndLinkers();
	virtual void CollapseLinkers();
	virtual bool CollapseNodes();
	virtual void SimplifyTopology();
	virtual void PostSimplification() {};
	int CheckCatenation();
	std::string GetCatenationInfo(int num_nets);

public:
	Deconstructor(OBMol* orig_mof);
	virtual ~Deconstructor() {};

	void SimplifyMOF(bool write_intermediate_cifs=true);

	// Output CIFs and building block identity.
	void SetOutputDir(const std::string &path);
	virtual void WriteCIFs();
	virtual std::string GetMOFInfo();

	// Utilities
	std::string GetOutputPath(const std::string &filename);
	void WriteSimplifiedNet(const std::string &base_filename);
	void WriteAtomsOfRole(const std::string &simplified_role, const std::string &base_filename="");
};


class MOFidDeconstructor : public Deconstructor {
// The original MOFid algorithm and implementation of Deconstructor (originally in sbu.cpp).
// Converts 4-c linkers in MIL-47, etc., to 2 x 3-c.
protected:
	virtual void PostSimplification();
	std::vector<std::string> PAsToUniqueInChIs(VirtualMol pa, const std::string &format);

public:
	MOFidDeconstructor(OBMol* orig_mof = NULL);
	virtual ~MOFidDeconstructor() {};
	std::string GetMOFkey(const std::string &topology = DEFAULT_MOFKEY_TOPOLOGY);
	std::string GetLinkerInChIs();
};


class SingleNodeDeconstructor : public Deconstructor {
// The single-node MOF deconstruction algorithm specified in 10.1021/acs.cgd.8b00126.
// Unlike the original deconstructor for MOFid building blocks, this algorithm detects the points
// of extension for the SBUs, which is more natural for tri-metallic clusters, etc., and always
// considers linkers as a single SBU (never any branch points).
protected:
	virtual void DetectInitialNodesAndLinkers();  // detect nodes as entire SBUs
	virtual void PostSimplification() {};  // do not inherit the MOFid 4-to-2x3 step
	void WriteSBUs(const std::string &base_filename, bool external_bond_pa, bool external_conn_pa);
	static std::pair<VirtualMol,VirtualMol> CalculateNonmetalRing(OBAtom* a, OBAtom* b);
	static VirtualMol GetNonmetalRingSubstituent(OBAtom* src);

public:
	SingleNodeDeconstructor(OBMol* orig_mof = NULL);
	virtual ~SingleNodeDeconstructor() {};
	virtual void WriteCIFs();
};


class AllNodeDeconstructor : public SingleNodeDeconstructor {
// The all-node MOF deconstruction algorithm specified in 10.1021/acs.cgd.8b00126.
// This algorithm is equivalent to the single node case, except for
// detecting branch points within linkers.
protected:
	// Branch and branch point identification in both the original MOF and simplified net pseudoatoms.
	// Need to use OBMol's for the pseudoatoms because the branches are converted to connectors,
	// and branch points may get simplified out of the net (e.g. solvent removal).
	VirtualMol branches_orig;
	VirtualMol branch_points_orig;
	OBMol branches_pa;
	OBMol branch_points_pa;

	virtual void CollapseLinkers();  // Detect branch points in the linker

	void TreeDecomposition(MappedMol *fragment_to_simplify, VirtualMol connection_points);
	void CollapseRings(MappedMol *fragment_to_simplify, bool fuse_bridged_rings = true);

public:
	AllNodeDeconstructor(OBMol* orig_mof = NULL);
	virtual ~AllNodeDeconstructor() {};
	virtual void WriteCIFs();
};


// Not implementing the ToposPro "standard representation" explicitly


} // end namespace OpenBabel
#endif // DECONSTRUCTOR_H

//! \file deconstructor.h
//! \brief deconstructor.h - Deconstruct a MOF into its building blocks and simplified net
