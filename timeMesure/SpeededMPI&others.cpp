#include <LinearSystemSolver.h>
#include <random>
#include <chrono>
#include <fstream>

#define TESTED_TYPE double

int main(){
    const int MATRIX_SIZE = 10;
    std::vector<TESTED_TYPE> a = {
        20,3,2,0,3,0,0,2,3,0,3,22,3,3,0,0,0,0,3,3,2,3,20,0,1,0,0,0,4,0,0,3,0,21,0,4,0,1,0,2,3,0,1,0,22,3,0,2,3,3,0,0,0,4,3,20,4,0,3,0,0,0,0,0,0,4,20,4,1,0,2,0,0,1,2,0,4,20,0,0,3,3,4,0,3,3,1,0,22,0,0,3,0,2,3,0,0,0,0,23

    };
    std::vector<TESTED_TYPE> x(MATRIX_SIZE);
    for(int i=0;i<MATRIX_SIZE;i++){
        x[i]=i;
    }
    const int N_TO_AV=5;
    
    CSRMatrix<TESTED_TYPE> A = CSRMatrix<TESTED_TYPE>::CSR_from_reg_matrix(Matrix<TESTED_TYPE>::create_matrix_from_array(a, MATRIX_SIZE, MATRIX_SIZE)); 
    std::cout<<"checkpoint\n";
    auto b = A*x;

    std::ofstream file;
    file.open("./TimeSpeededMPI.csv");
    file<<"N_iter, MPI_error, MPI_time, Speded_MPI_error, Speded_MPI_time, Yakobi_error, Yakobi_time, GaussSeidel_error, GaussSeidel_time\n";
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
            answ = ChebSpeededPIM(A, b, NIter, 1e-7, std::vector<TESTED_TYPE>(MATRIX_SIZE));
            for(int i=0; i<answ.size(); i++) {if(answ[i] - x[i] >= 1e-6) return 1;}
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