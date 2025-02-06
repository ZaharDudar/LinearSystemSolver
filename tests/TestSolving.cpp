#include <gtest/gtest.h>

#include <LinearSystemSolver.h>

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
