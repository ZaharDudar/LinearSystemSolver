#include <gtest/gtest.h>

#include <LinearSystemSolver.h>


TEST(CSRMatrixClass, CreateAndGetByIndx){
	std::map<std::pair<int, int>, int> doi={{{0,0},1},{{0,1},2},{{1,1},3},{{2,2},3}};
	CSRMatrix<int> tested_m(doi);

	EXPECT_EQ(tested_m(0,0), 1);
	EXPECT_EQ(tested_m(0,1), 2);
	EXPECT_EQ(tested_m(1,1), 3);
	EXPECT_EQ(tested_m(2,2), 3);
	EXPECT_EQ(tested_m(2,1), 0);
}

TEST(CSRMatrixClass, CSRFromRegular){
	std::vector<int> a{1,2,0,
					   4,5,6};
	Matrix<int> A = Matrix<int>::create_matrix_from_array(a, 2, 3);
	CSRMatrix<int> tested_m = CSRMatrix<int>::CSR_from_reg_matrix(A);

	EXPECT_EQ(tested_m(0,0), 1);
	EXPECT_EQ(tested_m(0,1), 2);
	EXPECT_EQ(tested_m(0,2), 0);
	EXPECT_EQ(tested_m(1,0), 4);
	EXPECT_EQ(tested_m(1,1), 5);
	EXPECT_EQ(tested_m(1,2), 6);

}

TEST(CSRMatrixClass, RegularFromCSR){
	std::map<std::pair<int, int>, int> doi={{{0,0},1},{{0,1},2},{{1,1},3},{{2,2},3}};
	CSRMatrix<int> tested_m(doi);

	Matrix<int> A = reg_matrix_from_CSR(tested_m);
	EXPECT_EQ(A(0,0), 1);
	EXPECT_EQ(A(0,1), 2);
	EXPECT_EQ(A(1,1), 3);
	EXPECT_EQ(A(2,2), 3);
	EXPECT_EQ(A(2,1), 0);
}

TEST(CSRMatrixClass, MultByVector1){
	std::vector<int> a{1,2,0,
					   4,5,6};
	Matrix<int> A = Matrix<int>::create_matrix_from_array(a, 2, 3);
	CSRMatrix<int> tested_m = CSRMatrix<int>::CSR_from_reg_matrix(A);

	std::vector<int> b{1,2,3};
	std::vector<int> true_answ{5,32};

	std::vector<int> getted_answ = tested_m * b;
	EXPECT_EQ(true_answ[0], getted_answ[0]);
	EXPECT_EQ(true_answ[1], getted_answ[1]);
	EXPECT_EQ(true_answ.size(), getted_answ.size());
}

TEST(CSRMatrixClass, MultByVector2){
	std::vector<int> a{1,2,0,
					   1,2,9,
					   7,2,0,
					   4,5,6};
	Matrix<int> A = Matrix<int>::create_matrix_from_array(a, 4, 3);
	CSRMatrix<int> tested_m = CSRMatrix<int>::CSR_from_reg_matrix(A);

	std::vector<int> b{1,2,3};
	std::vector<int> true_answ{5,32,11,32};

	std::vector<int> getted_answ = tested_m * b;
	EXPECT_EQ(true_answ[0], getted_answ[0]);
	EXPECT_EQ(true_answ[1], getted_answ[1]);
	EXPECT_EQ(true_answ[2], getted_answ[2]);
	EXPECT_EQ(true_answ[3], getted_answ[3]);
	EXPECT_EQ(true_answ.size(), getted_answ.size());
}