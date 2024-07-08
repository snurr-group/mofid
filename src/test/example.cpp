#include <vector>
#include "invector.h"
#include <gtest/gtest.h>

TEST(HelloTest, BasicAssertions) {
    // Expect two strings to not be equal
    EXPECT_STRNE("hello", "world");
    // Expect equality
    EXPECT_EQ(7 * 6, 42);
}

TEST(InVector, BasicAssertions) {
    std::vector<int> myVec{1, 2, 3};
    EXPECT_TRUE(inVector(1, myVec));
    EXPECT_TRUE(inVector(2, myVec));
    EXPECT_TRUE(inVector(3, myVec));
    EXPECT_FALSE(inVector(0, myVec));
    EXPECT_FALSE(inVector(4, myVec));
}
