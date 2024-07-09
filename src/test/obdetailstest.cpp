#include <gtest/gtest.h>
#include <array>
#include <openbabel/atom.h>

#include "obdetails.h"


TEST(IsMetalTest, HandlesMetal) {
}

TEST(IsMetalTest, HandlesNonMetal) {
    constexpr int numNonMetals{23};
    constexpr std::array<int, numNonMetals> nonMetalAtomicNums{1, 2, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18, 32, 33, 34, 35, 36, 52, 53, 54, 85, 86};
    OpenBabel::OBAtom nonMetalAtom{};
    for (int atomicNum : nonMetalAtomicNums) {
        nonMetalAtom.SetAtomicNum(atomicNum);
        const OpenBabel::OBAtom nonMetalAtom2{nonMetalAtom};
        const OpenBabel::OBAtom* nonMetalAtomPtr{&nonMetalAtom2};
        EXPECT_FALSE(OpenBabel::isMetal(nonMetalAtomPtr));
    }
}
