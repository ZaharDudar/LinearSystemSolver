#include <gtest/gtest.h>

#include <LinearSystemSolver.h>


TEST(VectorOps, Add){
    std::vector<int> a{1,2,3};
    std::vector<int> b{3,2,1};

    std::vector<int> c = a+b;
    EXPECT_EQ(c[0], 4); 
    EXPECT_EQ(c[1], 4);
    EXPECT_EQ(c[2], 4);
}

TEST(VectorOps, Sub){
    std::vector<int> a{1,2,3};
    std::vector<int> b{3,2,1};

    std::vector<int> c = a-b;
    EXPECT_EQ(c[0], -2); 
    EXPECT_EQ(c[1], 0);
    EXPECT_EQ(c[2], 2);
}

TEST(VectorOps, MulByValue){
    std::vector<int> a{1,2,3};

    std::vector<int> c = a*2;
    std::vector<int> d = 2*a;
    EXPECT_EQ(c[0], 2); 
    EXPECT_EQ(c[1], 4);
    EXPECT_EQ(c[2], 6);

    EXPECT_EQ(d[0], 2); 
    EXPECT_EQ(d[1], 4);
    EXPECT_EQ(d[2], 6);
}

TEST(VectorOps, Dot){
    std::vector<int> a{1,2,3};
    std::vector<int> b{3,2,1};

    int c = a*b;
    EXPECT_EQ(c, 10); 
}


