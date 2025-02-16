#include <gtest/gtest.h>

#include <LinearSystemSolver.h>

TEST(ThomasAlgorithm, Correct) {
	std::vector<float> a{1,2,0,
					     4,5,6,
					     0,8,9};

	std::vector<float> a_line{4,8};
	std::vector<float> b_line{1,5,9};
	std::vector<float> c_line{2,6};
	Matrix<float> A = Matrix<float>::create_matrix_from_array(a, 3, 3);
	std::vector<float> B{8,47,60};

	std::vector<float> right_solution{2, 3, 4};
	std::vector<float> getted_solution_dense = SolveThomasAglorithm(A,B);
	std::vector<float> getted_solution_threeDiag = SolveThomasAglorithm(a_line, b_line, c_line ,B);
	for(int i=0; i<3; i++){
		EXPECT_NEAR(getted_solution_dense[i],right_solution[i], 0.0001);
		EXPECT_NEAR(getted_solution_threeDiag[i],right_solution[i], 0.0001);
	}
	
	B =std::vector<float>{1,2,3};

	right_solution = std::vector<float>{0.04, 0.48, -0.0933};
    getted_solution_dense = SolveThomasAglorithm(A,B);
    getted_solution_threeDiag = SolveThomasAglorithm(a_line, b_line, c_line ,B);
	for(int i=0; i<3; i++){
		EXPECT_NEAR(getted_solution_dense[i],right_solution[i], 0.0001);
		EXPECT_NEAR(getted_solution_threeDiag[i],right_solution[i], 0.0001);
	}
}

TEST(QRHouseholder, test){
	std::vector<float> a{1,2,4,
						 3,3,2,
						 4,1,3};
	std::vector<float> b{1,2,3,
						4,5,6,
						4,5,6,
						7,8,9};
	std::vector<float> c{12, -51, 4,
						 6, 167,-68,
						 -4, 24, -41};
	std::vector<std::pair<std::pair<int, int>, std::vector<float>>> tests{{{3,3},a}, {{4,3},b}, {{3,3}, c}};

	for(auto &test : tests){
		Matrix<float> A = Matrix<float>::create_matrix_from_array(test.second,test.first.first, test.first.second);
		auto qr = hausholder_alorithm(A);
		auto recreated_A = qr.first * qr.second;
		for(int i=0; i<test.first.first; i++){
			for(int j=0; j<test.first.second; j++){
				EXPECT_NEAR(recreated_A(i,j), A(i,j), 0.0001);
			}
		}
	}

}

TEST(SolveWithQR, Correct){
	std::vector<float> a{1,2,4,
			3,3,2,
			4,1,3};
	std::vector<float> b{12, -51, 4,
			6, 167,-68,
			-4, 24, -41};
	std::vector<float> X{1,2,3};
	std::vector<std::pair<std::pair<int, int>, std::vector<float>>> tests{{{3,3},a}, {{3,3},b}};

	for(auto &test : tests){
		Matrix<float> A = Matrix<float>::create_matrix_from_array(test.second,test.first.first, test.first.second);
		auto qr = hausholder_alorithm(A);
		auto f = A * X;
		auto getted_answ = solve_by_QR(qr.first, qr.second, f);
		for(int i=0; i<3; i++){
			EXPECT_NEAR(X[i], getted_answ[i], 0.0001);
		}
	}
}
