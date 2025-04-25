#include <gtest/gtest.h>

#include <LinearSystemSolver.h>
// std::vector<float> a = {
//     4, 1, 2, 0, 3, 1,
//     1, 5, 0, 2, 1, 0,
//     2, 0, 6, 1, 0, 2,
//     0, 2, 1, 7, 1, 0,
//     3, 1, 0, 1, 8, 1,
//     1, 0, 2, 0, 1, 9
// };

TEST(GradientMethods, SteepestGradientDescent){
	std::vector<double> a{3,2,2,
		1,3,1,
		2,1,3};
	std::vector<double> X{1,2,3};
	std::vector<std::pair<std::pair<int, int>, std::vector<double>>> tests{{{3,3},a}};

	for(auto &test : tests){
		CSRMatrix<double> A = CSRMatrix<double>::CSR_from_reg_matrix(Matrix<double>::create_matrix_from_array(test.second, test.first.first, test.first.second));
		auto B = A * X;
		auto getted_answ = SteepestGradientDescent(A, B, 10000, (double)1e-7, std::vector<double>(3));
		for(int i=0; i<3; i++){
			EXPECT_NEAR(X[i], getted_answ[i], 0.0001);
		}
	}
}

TEST(GradientMethods, CG){
    std::vector<double> a = {
        4, 1, 2, 0, 3, 1,
        1, 5, 0, 2, 1, 0,
        2, 0, 6, 1, 0, 2,
        0, 2, 1, 7, 1, 0,
        3, 1, 0, 1, 8, 1,
        1, 0, 2, 0, 1, 9
    };
	std::vector<double> X{1,2,3,4,5,6};
	std::vector<std::pair<std::pair<int, int>, std::vector<double>>> tests{{{6,6},a}};

	for(auto &test : tests){
		CSRMatrix<double> A = CSRMatrix<double>::CSR_from_reg_matrix(Matrix<double>::create_matrix_from_array(test.second, test.first.first, test.first.second));
		auto B = A * X;
		auto getted_answ = CG(A, B, 10000, (double)1e-7, std::vector<double>(6));
		for(int i=0; i<6; i++){
			EXPECT_NEAR(X[i], getted_answ[i], 0.0001);
		}
	}
}
TEST(GradientMethods, BiCG){
	// std::vector<double> a{3,2,2,
	// 	1,3,1,
	// 	2,1,3};
	// std::vector<double> a = {
	// 	4, -1, -1, -0,  0,  0,
	// 	-1,  4, -0, -1,  0,  0,
	// 	-1, -0,  4, -1, -1, -0,
	// 	-0, -1, -1,  4, -0, -1,
	// 	0,  0, -1, -0,  4, -1,
	// 	0,  0, -0, -1, -1,  4,
	// };
	std::vector<double> a = {
        4, 1, 2, 0, 3, 1,
        1, 5, 0, 2, 1, 0,
        2, 0, 6, 1, 0, 2,
        0, 2, 1, 7, 1, 0,
        3, 1, 0, 1, 8, 1,
        1, 0, 2, 0, 1, 9
    };
	std::vector<double> X{1,2,3,4,5,6};
	std::vector<std::pair<std::pair<int, int>, std::vector<double>>> tests{{{6,6},a}};

	for(auto &test : tests){
		CSRMatrix<double> A = CSRMatrix<double>::CSR_from_reg_matrix(Matrix<double>::create_matrix_from_array(test.second, test.first.first, test.first.second));
		auto B = A * X;
		auto getted_answ = BiCG(A, B, 1000, (double)1e-7, std::vector<double>(6));
		for(int i=0; i<6; i++){
			EXPECT_NEAR(X[i], getted_answ[i], 0.0001);
		}
	}
}
TEST(GradientMethods, CGS){
    std::vector<double> a = {
        4, 1, 2, 0, 3, 1,
        1, 5, 0, 2, 1, 0,
        2, 0, 6, 1, 0, 2,
        0, 2, 1, 7, 1, 0,
        3, 1, 0, 1, 8, 1,
        1, 0, 2, 0, 1, 9
    };
	std::vector<double> X{1,2,3,4,5,6};
	std::vector<std::pair<std::pair<int, int>, std::vector<double>>> tests{{{6,6},a}};

	for(auto &test : tests){
		CSRMatrix<double> A = CSRMatrix<double>::CSR_from_reg_matrix(Matrix<double>::create_matrix_from_array(test.second, test.first.first, test.first.second));
		auto B = A * X;
		auto getted_answ = CGS(A, B, 100, (double)1e-7, std::vector<double>(6));
		for(int i=0; i<6; i++){
			EXPECT_NEAR(X[i], getted_answ[i], 0.0001);
		}
	}
}