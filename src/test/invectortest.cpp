#include <gtest/gtest.h>
#include <string>
#include <vector>

#include "invector.h"

using namespace std::string_literals; 

TEST(InVectorTest, HandlesIntVector) {
    std::vector<int> myVec{1, 2, 3};
    EXPECT_TRUE(inVector(1, myVec));
    EXPECT_TRUE(inVector(2, myVec));
    EXPECT_TRUE(inVector(3, myVec));
    EXPECT_FALSE(inVector(0, myVec));
    EXPECT_FALSE(inVector(4, myVec));
}

TEST(InVectorTest, HandlesBoolVector) {
    std::vector<bool> myVec{true, false};
    EXPECT_TRUE(inVector(true, myVec));
    EXPECT_TRUE(inVector(false, myVec));
    EXPECT_FALSE(inVector(true, {}));
    EXPECT_FALSE(inVector(false, {}));
}

TEST(InVectorTest, HandlesCharVector) {
    std::vector<char> myVec{'a', 'b', 'c'};
    EXPECT_TRUE(inVector('a', myVec));
    EXPECT_TRUE(inVector('b', myVec));
    EXPECT_TRUE(inVector('c', myVec));
    EXPECT_FALSE(inVector('d', myVec));
    EXPECT_FALSE(inVector('e', myVec));
}

TEST(InVectorTest, HandlesDoubleVector) {
    std::vector<double> myVec{1.1, 2.2, 3.3};
    EXPECT_TRUE(inVector(1.1, myVec));
    EXPECT_TRUE(inVector(2.2, myVec));
    EXPECT_TRUE(inVector(3.3, myVec));
    EXPECT_FALSE(inVector(0.0, myVec));
    EXPECT_FALSE(inVector(4.4, myVec));
}

TEST(InVectorTest, HandlesStringVector) {
    std::vector<std::string> myVec{"apple", "banana", "cherry"};
    EXPECT_TRUE(inVector("apple"s, myVec));
    EXPECT_TRUE(inVector("banana"s, myVec));
    EXPECT_TRUE(inVector("cherry"s, myVec));
    EXPECT_FALSE(inVector("date"s, myVec));
    EXPECT_FALSE(inVector("fig"s, myVec));
}
