#include <gtest/gtest.h>

#include <LinearSystemSolver.h>
// Demonstrate some basic assertions.
TEST(MatrixClass, DiagGen) {
	std::vector<int> diag{1,2,3};
	Matrix<int> test_mtx = Matrix<int>::create_diag_matrix(diag);
	EXPECT_EQ(test_mtx.shape(), std::make_pair(3,3));


	for(int i=0; i<3; i++){
		for(int j=0; j<3 ; j++){
			if(i==j){
				EXPECT_EQ(test_mtx(i,j), diag[i]);
			}
			else{
				EXPECT_EQ(test_mtx(i,j), 0);
			}
		}
	}
}

TEST(MatrixClass, ThreeDiagGen) {
	std::vector<int> diagA{1,1};
	std::vector<int> diagB{2,2,2};
	std::vector<int> diagC{3,3};
	EXPECT_ANY_THROW(Matrix<int>::create_3diag_matrix(diagA, diagC, diagA)) << "Does'n throw";
	Matrix<int> test_mtx = Matrix<int>::create_3diag_matrix(diagA, diagB, diagC);
	EXPECT_EQ(test_mtx.shape(), std::make_pair(3,3)) << "Wrong shape!";

    for(int i=0; i<2; i++){
        EXPECT_EQ(test_mtx(i,i+1), diagA[i]) << "Wrong diagA";
    }
    for(int i=0; i<2; i++){
		EXPECT_EQ(test_mtx(i+1,i), diagC[i]) << "Wrong diagC";
    }
}

