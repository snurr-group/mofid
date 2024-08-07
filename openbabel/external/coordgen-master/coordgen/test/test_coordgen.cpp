#define BOOST_TEST_MODULE Test_Coordgen

#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>
#include <unordered_set>

#include "../sketcherMinimizer.h"
#include "../sketcherMinimizerMaths.h"
#include "../sketcherMaeReading.h"

#include "maeparser/MaeConstants.hpp"
#include "maeparser/Reader.hpp"

using std::unordered_set;
using namespace schrodinger;

const boost::filesystem::path test_samples_path(TEST_SAMPLES_PATH);

namespace
{
std::map<sketcherMinimizerAtom*, int>
getReportingIndices(sketcherMinimizerMolecule& mol)
{
    std::map<sketcherMinimizerAtom*, int> fakeIndices;
    int index = 0;
    for (auto& atom : mol.getAtoms()) {
        fakeIndices.emplace(atom, ++index);
    }
    return fakeIndices;
}

bool areBondsNearIdeal(sketcherMinimizerMolecule& mol,
                       std::map<sketcherMinimizerAtom*, int>& indices)
{
    const float targetBondLength = BONDLENGTH * BONDLENGTH;
    const auto tolerance = static_cast<float>(targetBondLength * 0.1);

    bool passed = true;
    for (auto& bond : mol.getBonds()) {
        auto& startCoordinates = bond->getStartAtom()->getCoordinates();
        auto& endCoordinates = bond->getEndAtom()->getCoordinates();

        const auto sq_distance = sketcherMinimizerMaths::squaredDistance(
            startCoordinates, endCoordinates);
        const auto deviation = sq_distance - targetBondLength;
        if (deviation < -tolerance || deviation > tolerance) {
            std::cerr << "Bond" << indices[bond->getStartAtom()] << '-'
                      << indices[bond->getEndAtom()] << " has length "
                      << sq_distance << " (" << targetBondLength << ")\n";
            passed = false;
        }
    }

    return passed;
}
} // namespace

BOOST_AUTO_TEST_CASE(SampleTest)
{
    // A small sample test showing how to import a molecule from a .mae file and
    // initialize the minimizer with it.

    const std::string testfile = (test_samples_path / "test.mae").string();

    mae::Reader r(testfile);
    auto block = r.next(mae::CT_BLOCK);
    BOOST_REQUIRE(block != nullptr);

    auto* mol = mol_from_mae_block(*block);
    BOOST_REQUIRE(mol != nullptr);
    BOOST_REQUIRE_EQUAL(mol->getAtoms().size(), 26);
    BOOST_REQUIRE_EQUAL(mol->getBonds().size(), 26);

    sketcherMinimizer minimizer;
    minimizer.initialize(mol); // minimizer takes ownership of mol
    minimizer.runGenerateCoordinates();

    for (auto& atom : mol->getAtoms()) {
        auto c = atom->getCoordinates();

        // This is best we can do with checking: there are precision and
        // rounding issues depending on platform and environment.
        BOOST_CHECK(c.x() != 0 || c.y() != 0);
    }

    auto indices = getReportingIndices(*mol);
    BOOST_CHECK(areBondsNearIdeal(*mol, indices));
}


BOOST_AUTO_TEST_CASE(TemplateTest)
{
    ///
    // Do the structures in the templates file get the same coordinates that
    // were supplied in the templates file?

    const boost::filesystem::path source_dir(SOURCE_DIR);
    const std::string templates_file = (source_dir / "templates.mae").string();

    // Known issues. See issue #52
    const unordered_set<size_t> no_match = {1, 8, 19, 20, 22, 32, 43, 53, 65, 66, 67};
    // 32 is odd. minimization removes atoms?? But it matches??
    const unordered_set<size_t> match_incorrectly = {18, 27};

    mae::Reader r(templates_file);
    std::shared_ptr<mae::Block> b;
    size_t template_index = 0;
    while ((b = r.next(mae::CT_BLOCK)) != nullptr) {
        if (no_match.count(template_index) > 0) {
            ++template_index;
            continue;
        }

        auto* mol = mol_from_mae_block(*b);
        BOOST_REQUIRE(mol != nullptr);
        const auto original_atom_count = mol->getAtoms().size();

        sketcherMinimizer minimizer;
        minimizer.setTemplateFileDir(source_dir.string());

        minimizer.initialize(mol); // minimizer takes ownership of mol
        minimizer.runGenerateCoordinates();

        BOOST_CHECK_EQUAL(original_atom_count, mol->getAtoms().size());

        bool any_rigid = false;
        bool all_rigid = true;
        for (auto a: mol->getAtoms()) {
            if (a->rigid) {
                any_rigid = true;
            } else {
                all_rigid = false;
            }
        }
        const bool matches_incorrectly = match_incorrectly.count(template_index) > 0;
        if (!any_rigid) {
            BOOST_CHECK_MESSAGE(any_rigid, "No template found for " << template_index);
        } else if (!matches_incorrectly) {
            BOOST_CHECK_MESSAGE(all_rigid, "Not all atoms templated for " << template_index);
        }

        ++template_index;
    }
}

BOOST_AUTO_TEST_CASE(ClearWedgesTest)
{
    // test that when writing stereochemistry we first clear the existing one

    const std::string testfile = (test_samples_path / "testChirality.mae").string();

    mae::Reader r(testfile);
    auto block = r.next(mae::CT_BLOCK);
    BOOST_REQUIRE(block != nullptr);

    auto* mol = mol_from_mae_block(*block);
    BOOST_REQUIRE(mol != nullptr);

     /*set wedges on all bonds*/
    mol->getBonds().at(0)->hasStereochemistryDisplay = true;
    mol->getBonds().at(1)->hasStereochemistryDisplay = true;
    mol->getBonds().at(2)->hasStereochemistryDisplay = true;
    mol->getBonds().at(3)->hasStereochemistryDisplay = true;

    /*set chirality on the atom*/
    auto carbon = mol->getAtoms().at(0);
    BOOST_REQUIRE_EQUAL(carbon->atomicNumber, 6);
    carbon->hasStereochemistrySet = true;
    carbon->isR = true;

    /*run minimization*/
    sketcherMinimizer minimizer;
    minimizer.initialize(mol); // minimizer takes ownership of mol
    minimizer.runGenerateCoordinates();

    /*make sure that previous wedges are reset and only 2 are now used*/
    BOOST_REQUIRE_EQUAL(mol->getBonds().at(0)->hasStereochemistryDisplay, false);
    BOOST_REQUIRE_EQUAL(mol->getBonds().at(1)->hasStereochemistryDisplay, true);
    BOOST_REQUIRE_EQUAL(mol->getBonds().at(2)->hasStereochemistryDisplay, false);
    BOOST_REQUIRE_EQUAL(mol->getBonds().at(3)->hasStereochemistryDisplay, true);
}
