#include <gtest/gtest.h>

#include <LinearSystemSolver.h>
#include <map>

//Check creating diagonal matrix from vector
TEST(ColMajMatrixClass, DiagGen) {
	std::vector<int> diag{1,2,3};
	ColMajMatrix<int> test_mtx = ColMajMatrix<int>::create_diag_matrix(diag);
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
TEST(ColMajMatrixClass, ThreeDiagGen) {
	std::vector<int> diagA{1,1};
	std::vector<int> diagB{2,2,2};
	std::vector<int> diagC{3,3};
	EXPECT_ANY_THROW(ColMajMatrix<int>::create_3diag_matrix(diagA, diagC, diagA)) << "Doesn't throw";
	ColMajMatrix<int> test_mtx = ColMajMatrix<int>::create_3diag_matrix(diagA, diagB, diagC);
	EXPECT_EQ(test_mtx.shape(), std::make_pair(3,3)) << "Wrong shape!";

    for(int i=0; i<2; i++){
        EXPECT_EQ(test_mtx(i+1,i), diagA[i]) << "Wrong diagA";
    }
    for(int i=0; i<2; i++){
		EXPECT_EQ(test_mtx(i,i+1), diagC[i]) << "Wrong diagC";
    }
}

//Check for creating matrix from one-dim array with given shape
TEST(ColMajMatrixClass, MatrixFromArray) {
	std::vector<int> a{1,2,3,
					   4,5,6,
					   7,8,9};
	EXPECT_ANY_THROW(ColMajMatrix<int>::create_matrix_from_array(a, 2, 3)) << "Doesn't throw";
	ColMajMatrix<int> test_mtx = ColMajMatrix<int>::create_matrix_from_array(a, 3, 3);

	EXPECT_EQ(test_mtx.shape(), std::make_pair(3,3)) << "Wrong shape!";

	for(int i=0; i<3; i++){
		for(int j=0; j<3 ; j++){
			EXPECT_EQ(test_mtx(i,j), a[i*3+j]);
		}
	}
	EXPECT_EQ(test_mtx(1,0), 4);


	a = {1,2,0,
		 4,5,6};
	EXPECT_ANY_THROW(ColMajMatrix<int>::create_matrix_from_array(a, 3, 3)) << "Doesn't throw";
	test_mtx = ColMajMatrix<int>::create_matrix_from_array(a, 2, 3);
	EXPECT_EQ(test_mtx(0,0), 1);
	EXPECT_EQ(test_mtx(0,1), 2);
	EXPECT_EQ(test_mtx(0,2), 0);
	EXPECT_EQ(test_mtx(1,0), 4);
	EXPECT_EQ(test_mtx(1,1), 5);
	EXPECT_EQ(test_mtx(1,2), 6);
}

TEST(ColMajMatrixClass, MultByVector1){
	std::vector<int> a{1,2,0,
					   4,5,6};
	ColMajMatrix<int> A = ColMajMatrix<int>::create_matrix_from_array(a, 2, 3);

	std::vector<int> b{1,2,3};
	std::vector<int> true_answ{5,32};

	std::vector<int> getted_answ = A * b;
	EXPECT_EQ(true_answ[0], getted_answ[0]);
	EXPECT_EQ(true_answ[1], getted_answ[1]);
	EXPECT_EQ(true_answ.size(), getted_answ.size());
}

TEST(ColMajMatrixClass, MultByVector2){
	std::vector<int> a{1,2,0,
					   1,2,9,
					   7,2,0,
					   4,5,6};
	ColMajMatrix<int> A = ColMajMatrix<int>::create_matrix_from_array(a, 4, 3);

	std::vector<int> b{1,2,3};
	std::vector<int> true_answ{5,32,11,32};

	std::vector<int> getted_answ = A * b;
	EXPECT_EQ(true_answ[0], getted_answ[0]);
	EXPECT_EQ(true_answ[1], getted_answ[1]);
	EXPECT_EQ(true_answ[2], getted_answ[2]);
	EXPECT_EQ(true_answ[3], getted_answ[3]);
	EXPECT_EQ(true_answ.size(), getted_answ.size());
}
