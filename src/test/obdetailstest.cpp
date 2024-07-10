#include <gtest/gtest.h>
#include <array>
#include <map>

#include "obdetails.h"

#include <openbabel/atom.h>

namespace OBDetailsTest {
    constexpr int numElements{118};
    constexpr int numNonMetals{23};
    constexpr std::array<int, numNonMetals> nonMetalAtomicNums{1, 2, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18, 32, 33, 34, 35, 36, 52, 53, 54, 85, 86};
    const std::map<int, const char*> elements{{1, "Hydro"}, {2, "Heliu"}, {3, "Lithi"}, {8, "Oxyge"}, {14, "Silic"}, {29, "Coppe"}, {39, "Yttri"}, {48, "Cadmi"}, {54, "Xenon"}, {81, "Thall"}, {101, "Mende"}};
}

using namespace OBDetailsTest;

TEST(IsMetalTest, HandlesMetal) {
    int idx{0};
    OpenBabel::OBAtom metalAtom{};
    for (std::size_t atomicNum{1}; atomicNum <= numElements; ++atomicNum) {
        if ((idx < numNonMetals) && (atomicNum == nonMetalAtomicNums[idx]))
            ++idx;
        else {
            metalAtom.SetAtomicNum(atomicNum);
            const OpenBabel::OBAtom metalAtom2{metalAtom};
            EXPECT_TRUE(OpenBabel::isMetal(&metalAtom2));
        }
    }
}

TEST(IsMetalTest, HandlesNonMetal) {
    OpenBabel::OBAtom nonMetalAtom{};
    for (int atomicNum : nonMetalAtomicNums) {
        nonMetalAtom.SetAtomicNum(atomicNum);
        const OpenBabel::OBAtom nonMetalAtom2{nonMetalAtom};
        EXPECT_FALSE(OpenBabel::isMetal(&nonMetalAtom2));
    }
}

TEST (ChangeAtomElementTest, HandlesAtom) {
    OpenBabel::OBAtom atom{};
    for (const auto& p : OBDetailsTest::elements) {
        const int atomicNum{p.first};
        const char* type{p.second};
        OpenBabel::changeAtomElement(&atom, atomicNum);
        EXPECT_EQ(atomicNum, atom.GetAtomicNum());
        EXPECT_STREQ(type, atom.GetType());
    }
}

TEST (AtomsEqualTest, HandlesUnequal) {
}

TEST (AtomsEqualTest, HandlesEqual) {
}

TEST (RTrimWhiteSpaceTest, HandlesEmptyString) {
    EXPECT_EQ("", OpenBabel::rtrimWhiteSpace(""));
}

TEST (RTrimWhiteSpaceTest, HandlesTrimmedString) {
    EXPECT_EQ("Hi", OpenBabel::rtrimWhiteSpace("Hi"));
    EXPECT_EQ("Bye", OpenBabel::rtrimWhiteSpace("Bye"));
}

TEST (RTrimWhiteSpaceTest, HandlesUntrimmedString) {
    EXPECT_EQ("Hi", OpenBabel::rtrimWhiteSpace("Hi       "));
    EXPECT_EQ("Bye", OpenBabel::rtrimWhiteSpace("Bye "));
    EXPECT_EQ("     Hi", OpenBabel::rtrimWhiteSpace("     Hi"));
    EXPECT_EQ(" Bye", OpenBabel::rtrimWhiteSpace(" Bye"));
}