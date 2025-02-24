#include <LinearSystemSolver.h>
#include <random>
#include <chrono>
#include <fstream>

#define TESTED_TYPE double

int main(){
    const int MATRIX_SIZE = 6;
    std::vector<TESTED_TYPE> a = {
        4, 1, 2, 0, 3, 1,
        1, 5, 0, 2, 1, 0,
        2, 0, 6, 1, 0, 2,
        0, 2, 1, 7, 1, 0,
        3, 1, 0, 1, 8, 1,
        1, 0, 2, 0, 1, 9
    };
    std::vector<TESTED_TYPE> x{1,2,3,4,5,6};
    const int N_TO_AV=50;

    CSRMatrix<TESTED_TYPE> A = CSRMatrix<TESTED_TYPE>::CSR_from_reg_matrix(Matrix<TESTED_TYPE>::create_matrix_from_array(a, 6,6)); 
    auto b = A*x;

    std::ofstream file;
    file.open("./TimeIterMethods.csv");
    file<<"N_iter, MPI_error, MPI_time, Yakobi_error, Yakobi_time, GaussSeidel_error, GaussSeidel_time\n";
    auto st = std::chrono::high_resolution_clock::now();
    auto en = std::chrono::high_resolution_clock::now();
    double time;
    std::vector<TESTED_TYPE> answ;
    double error;
    for(int NIter = 1; NIter <= 1000; NIter+=1){
        file<<NIter<<",";

        error = 0;
        time = 0;
        for(int a=0; a < N_TO_AV; a++){
            st = std::chrono::high_resolution_clock::now();
            answ = PrimeIterMethod(A, b, NIter, 1e-7, std::vector<TESTED_TYPE>(MATRIX_SIZE));
            en = std::chrono::high_resolution_clock::now();
            error += abs(answ - x);
            time += std::chrono::duration_cast<std::chrono::microseconds>(en-st).count();
        }
        file<<error/N_TO_AV<<","<<time/N_TO_AV<<",";

        error = 0;
        time = 0;
        for(int a=0; a < N_TO_AV; a++){
            st = std::chrono::high_resolution_clock::now();
            answ = YakobyMethod(A, b, NIter, 1e-7, std::vector<TESTED_TYPE>(MATRIX_SIZE));
            en = std::chrono::high_resolution_clock::now();
            error += abs(answ - x);
            time += std::chrono::duration_cast<std::chrono::microseconds>(en-st).count();
        }
        file<<error/N_TO_AV<<","<<time/N_TO_AV<<",";

        error = 0;
        time = 0;
        for(int a=0; a < N_TO_AV; a++){
            st = std::chrono::high_resolution_clock::now();
            answ = GaussSeidelMethod(A, b, NIter, 1e-7, std::vector<TESTED_TYPE>(MATRIX_SIZE));
            en = std::chrono::high_resolution_clock::now();
            error += abs(answ - x);
            time += std::chrono::duration_cast<std::chrono::microseconds>(en-st).count();
        }
        file<<error/N_TO_AV<<","<<time/N_TO_AV<<"\n";        
    }
    file.close();

}