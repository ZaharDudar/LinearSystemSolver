#include <gtest/gtest.h>

#include <LinearSystemSolver.h>
#include <IterMethods.h>

TEST(IterMethods, PrimeIter_Method){
	std::vector<double> a{3,2,2,
			1,3,1,
			2,1,3};
	std::vector<double> X{1,2,3};
	std::vector<std::pair<std::pair<int, int>, std::vector<double>>> tests{{{3,3},a}};

	for(auto &test : tests){
		CSRMatrix<double> A = CSRMatrix<double>::CSR_from_reg_matrix(Matrix<double>::create_matrix_from_array(test.second, test.first.first, test.first.second));
		auto B = A * X;
		auto getted_answ = PrimeIterMethod(A, B, 10000, (double)1e-7, std::vector<double>(3));
		for(int i=0; i<3; i++){
			EXPECT_NEAR(X[i], getted_answ[i], 0.0001);
		}
	}
}

TEST(IterMethods, Yakoby_Method){
	std::vector<float> a{9,2,2,
        1,9,1,
        2,1,9};
    std::vector<float> X{1,2,3};
	std::vector<std::pair<std::pair<int, int>, std::vector<float>>> tests{{{3,3},a}};

	for(auto &test : tests){
		CSRMatrix<float> A = CSRMatrix<float>::CSR_from_reg_matrix(Matrix<float>::create_matrix_from_array(test.second, test.first.first, test.first.second));
		auto B = A * X;
		auto getted_answ = JacobiMethod(A, B, 1000, (float)1e-7, std::vector<float>(3));
		for(int i=0; i<3; i++){
			EXPECT_NEAR(X[i], getted_answ[i], 0.0001);
		}
	}
}

TEST(IterMethods, GaussSeidel_Method){
    std::vector<float> a = {
        4, 1, 2, 0, 3, 1,
        1, 5, 0, 2, 1, 0,
        2, 0, 6, 1, 0, 2,
        0, 2, 1, 7, 1, 0,
        3, 1, 0, 1, 8, 1,
        1, 0, 2, 0, 1, 9
    };

	std::vector<float> X{1,2,3,4,5,6};
	std::vector<std::pair<std::pair<int, int>, std::vector<float>>> tests{{{6,6},a}};

	for(auto &test : tests){
		CSRMatrix<float> A = CSRMatrix<float>::CSR_from_reg_matrix(Matrix<float>::create_matrix_from_array(test.second, test.first.first, test.first.second));
		auto B = A * X;
		auto getted_answ = GaussSeidelMethod(A, B, 1000, (float)1e-7, std::vector<float>(6));
		for(int i=0; i<6; i++){
			EXPECT_NEAR(X[i], getted_answ[i], 0.0001);
		}
	}
}

TEST(IterMethods, SpeededPrimeIter_Method){
	std::vector<double> a{3,2,2,
						  1,3,1,
						  2,1,3};
	std::vector<double> X{1,2,3};
	std::vector<std::pair<std::pair<int, int>, std::vector<double>>> tests{{{3,3},a}};

	for(auto &test : tests){
		CSRMatrix<double> A = CSRMatrix<double>::CSR_from_reg_matrix(Matrix<double>::create_matrix_from_array(test.second, test.first.first, test.first.second));
		auto B = A * X;
		auto getted_answ = ChebSpeededPIM(A, B, 10000, (double)1e-7, std::vector<double>(3));
		for(int i=0; i<3; i++){
			EXPECT_NEAR(X[i], getted_answ[i], 0.0001);
		}
	}
}



TEST(IterMethods, SuccessiveOverRelaxation){
    std::vector<float> a = {
        4, 1, 2, 0, 3, 1,
        1, 5, 0, 2, 1, 0,
        2, 0, 6, 1, 0, 2,
        0, 2, 1, 7, 1, 0,
        3, 1, 0, 1, 8, 1,
        1, 0, 2, 0, 1, 9
    };

	std::vector<float> X{1,2,3,4,5,6};
	std::vector<std::pair<std::pair<int, int>, std::vector<float>>> tests{{{6,6},a}};

	for(auto &test : tests){
		CSRMatrix<float> A = CSRMatrix<float>::CSR_from_reg_matrix(Matrix<float>::create_matrix_from_array(test.second, test.first.first, test.first.second));
		auto B = A * X;
		auto getted_answ = SOR(A, B, 1000, (float)1e-7, std::vector<float>(6));
		for(int i=0; i<6; i++){
			EXPECT_NEAR(X[i], getted_answ[i], 0.0001);
		}
	}
}

TEST(IterMethods, SymmetricGSMethod){
    std::vector<float> a = {
        4, 1, 2, 0, 3, 1,
        1, 5, 0, 2, 1, 0,
        2, 0, 6, 1, 0, 2,
        0, 2, 1, 7, 1, 0,
        3, 1, 0, 1, 8, 1,
        1, 0, 2, 0, 1, 9
    };

	std::vector<float> X{1,2,3,4,5,6};
	std::vector<std::pair<std::pair<int, int>, std::vector<float>>> tests{{{6,6},a}};

	for(auto &test : tests){
		CSRMatrix<float> A = CSRMatrix<float>::CSR_from_reg_matrix(Matrix<float>::create_matrix_from_array(test.second, test.first.first, test.first.second));
		auto B = A * X;
		auto getted_answ = SymmetricGSMethod(A, B, 1000, (float)1e-7, std::vector<float>(6));
		for(int i=0; i<6; i++){
			EXPECT_NEAR(X[i], getted_answ[i], 0.0001);
		}
	}
}

TEST(UniversalBoost, PrimeIter_Method){
	std::vector<double> a{	3,1,2,
							1,3,1,
							2,1,3};
	std::vector<double> X{1,2,3};
	std::vector<std::pair<std::pair<int, int>, std::vector<double>>> tests{{{3,3},a}};

	for(auto &test : tests){
		CSRMatrix<double> A = CSRMatrix<double>::CSR_from_reg_matrix(Matrix<double>::create_matrix_from_array(test.second, test.first.first, test.first.second));
		auto B = A * X;
		auto step = [](CSRMatrix<double> &A, std::vector<double> &b, std::vector<double> &x0) { return PrimeIterMethod(A, b, 1, 1e-7, x0);};
		auto getted_answ = universalChebBoost<double>(step, A, B, 1000, (double)1e-7,std::vector<double>(3));
		for(int i=0; i<3; i++){
			EXPECT_NEAR(X[i], getted_answ[i], 0.0001);
		}
	}
}

TEST(UniversalBoost, Jacobi_Method){
	std::vector<float> a{9,2,2,
        1,9,1,
        2,1,9};
    std::vector<float> X{1,2,3};
	std::vector<std::pair<std::pair<int, int>, std::vector<float>>> tests{{{3,3},a}};

	for(auto &test : tests){
		CSRMatrix<float> A = CSRMatrix<float>::CSR_from_reg_matrix(Matrix<float>::create_matrix_from_array(test.second, test.first.first, test.first.second));
		auto B = A * X;
		auto step = [](CSRMatrix<float> &A, std::vector<float> &b, std::vector<float> &x0) { return JacobiMethod(A, b, 1, (float)1e-7, x0);};
		auto getted_answ = universalChebBoost<float>(step, A, B, 1000, (float)1e-7,std::vector<float>(3));
		for(int i=0; i<3; i++){
			EXPECT_NEAR(X[i], getted_answ[i], 0.0001);
		}
	}
}

TEST(UniversalBoost, GaussSeidel_Method){
    // std::vector<float> a = {
    //     4, 1, 2, 0, 3, 1,
    //     1, 5, 0, 2, 1, 0,
    //     2, 0, 6, 1, 0, 2,
    //     0, 2, 1, 7, 1, 0,
    //     3, 1, 0, 1, 8, 1,
    //     1, 0, 2, 0, 1, 9
    // };
    std::vector<float> a = {
		 4, -1, -1, -0,  0,  0,
		-1,  4, -0, -1,  0,  0,
		-1, -0,  4, -1, -1, -0,
		-0, -1, -1,  4, -0, -1,
		 0,  0, -1, -0,  4, -1,
		 0,  0, -0, -1, -1,  4,
    };

	std::vector<float> X{1,2,3,4,5,6};
	std::vector<std::pair<std::pair<int, int>, std::vector<float>>> tests{{{6,6},a}};

	for(auto &test : tests){
		CSRMatrix<float> A = CSRMatrix<float>::CSR_from_reg_matrix(Matrix<float>::create_matrix_from_array(test.second, test.first.first, test.first.second));
		auto B = A * X;
		auto step = [](CSRMatrix<float> &A, std::vector<float> &b, std::vector<float> &x0) { return GaussSeidelMethod(A, b, 1, (float)1e-7, x0);};
		auto getted_answ = universalChebBoost<float>(step, A, B, 1000, (float)1e-7, std::vector<float>(6));

		for(int i=0; i<6; i++){
			EXPECT_NEAR(X[i], getted_answ[i], 0.0001);
		}
	}
}