#include <gtest/gtest.h>

#include <LinearSystemSolver.h>
//Check creating diagonal matrix from vector
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

//Check creating array with main and 2 nearest diags
TEST(MatrixClass, ThreeDiagGen) {
	std::vector<int> diagA{1,1};
	std::vector<int> diagB{2,2,2};
	std::vector<int> diagC{3,3};
	EXPECT_ANY_THROW(Matrix<int>::create_3diag_matrix(diagA, diagC, diagA)) << "Doesn't throw";
	Matrix<int> test_mtx = Matrix<int>::create_3diag_matrix(diagA, diagB, diagC);
	EXPECT_EQ(test_mtx.shape(), std::make_pair(3,3)) << "Wrong shape!";

    for(int i=0; i<2; i++){
        EXPECT_EQ(test_mtx(i+1,i), diagA[i]) << "Wrong diagA";
    }
    for(int i=0; i<2; i++){
		EXPECT_EQ(test_mtx(i,i+1), diagC[i]) << "Wrong diagC";
    }
}

//Check for creating matrix from one-dim array with given shape
TEST(MatrixClass, MatrixFromArray) {
	std::vector<int> a{1,2,3,
					   4,5,6,
					   7,8,9};
	EXPECT_ANY_THROW(Matrix<int>::create_matrix_from_array(a, 2, 3)) << "Doesn't throw";
	Matrix<int> test_mtx = Matrix<int>::create_matrix_from_array(a, 3, 3);
	EXPECT_EQ(test_mtx.shape(), std::make_pair(3,3)) << "Wrong shape!";

	for(int i=0; i<3; i++){
		for(int j=0; j<3 ; j++){
			EXPECT_EQ(test_mtx(i,j), a[i*3+j]);
		}
	}
	EXPECT_EQ(test_mtx(1,0), 4);
}


TEST(ThomasAlgorithm, Correct) {
	std::vector<float> a{1,2,0,
					     4,5,6,
					     0,8,9};
	Matrix<float> A = Matrix<float>::create_matrix_from_array(a, 3, 3);
	std::vector<float> B{8,47,60};

	std::vector<float> right_solution{2, 3, 4};
	std::vector<float> getted_solution = SolveThomasAglorithm(A,B);
	for(int i=0; i<3; i++){
		EXPECT_NEAR(getted_solution[i],right_solution[i], 0.0001);
	}
	
	B =std::vector<float>{1,2,3};

	right_solution = std::vector<float>{0.04, 0.48, -0.0933};
    getted_solution = SolveThomasAglorithm(A,B);
	for(int i=0; i<3; i++){
		EXPECT_NEAR(getted_solution[i],right_solution[i], 0.0001);
	}
}

