#include <LinearSystemSolver.h>

int main(){
    std::vector<int> diagA{1,1};
	std::vector<int> diagB{2,2,2};
	std::vector<int> diagC{3,3};
	Matrix<int>::create_3diag_matrix(diagA, diagC, diagA);
    return 0;
}