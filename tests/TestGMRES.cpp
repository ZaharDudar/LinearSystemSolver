#include <gtest/gtest.h>

#include <LinearSystemSolver.h>

TEST(GMRES, FirstTest){
	std::vector<double> a{1,2,0,
						 2,1,3,
						 0,0,1};
    std::vector<double> X{1,2,3};
	std::vector<std::pair<std::pair<int, int>, std::vector<double>>> tests{{{3,3},a}};

	for(auto &test : tests){
		CSRMatrix<double> A = CSRMatrix<double>::CSR_from_reg_matrix(Matrix<double>::create_matrix_from_array(test.second, test.first.first, test.first.second));
		auto B = A * X;
		auto getted_answ = GMRES(A, B, 1000, 20, (double)1e-7, std::vector<double>(3));
		for(int i=0; i<3; i++){
			EXPECT_NEAR(X[i], getted_answ[i], 0.0001);
		}
	}
}