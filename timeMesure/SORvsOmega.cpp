#include <LinearSystemSolver.h>
#include <random>
#include <chrono>
#include <fstream>

#define TESTED_TYPE double

TESTED_TYPE get_omega(int omega){return omega*0.1+1;}

int main(){
    const int MATRIX_SIZE = 10;
    std::vector<TESTED_TYPE> a = {
        20,3,2,0,3,0,0,2,3,0,3,22,3,3,0,0,0,0,3,3,2,3,20,0,1,0,0,0,4,0,0,3,0,21,0,4,0,1,0,2,3,0,1,0,22,3,0,2,3,3,0,0,0,4,3,20,4,0,3,0,0,0,0,0,0,4,20,4,1,0,2,0,0,1,2,0,4,20,0,0,3,3,4,0,3,3,1,0,22,0,0,3,0,2,3,0,0,0,0,23

    };
    std::vector<TESTED_TYPE> x(MATRIX_SIZE);
    for(int i=0;i<MATRIX_SIZE;i++){
        x[i]=i;
    }
    const int N_TO_AV=15;
    const int MAX_OMEGA_INT=10;


    CSRMatrix<TESTED_TYPE> A = CSRMatrix<TESTED_TYPE>::CSR_from_reg_matrix(Matrix<TESTED_TYPE>::create_matrix_from_array(a, MATRIX_SIZE, MATRIX_SIZE)); 
    auto b = A*x;
    
    std::ofstream file;
    file.open("./SORvsOmega.csv");
    file<<"N_iter,";
    std::cout<<"checkpoint\n";
    for(int omega = 0; omega<=MAX_OMEGA_INT; omega++){
        file<<"SOR_"<<get_omega(omega)<<"_error,";
        file<<"SOR_"<<get_omega(omega)<<"_time,";
    }
    file<<"SOR_optimal_error,";
    file<<"SOR_optimal_time,";
    file<<"\n";
    auto st = std::chrono::high_resolution_clock::now();
    auto en = std::chrono::high_resolution_clock::now();
    double time;
    std::vector<TESTED_TYPE> answ;
    double error;
    for(int NIter = 1; NIter <= 1000; NIter+=1){
        file<<NIter<<",";
        for(int omega = 0; omega<=MAX_OMEGA_INT; omega++){

            error = 0;
            time = 0;
            for(int a=0; a < N_TO_AV; a++){
                st = std::chrono::high_resolution_clock::now();
                answ = SOR(A, b, NIter, 1e-7, std::vector<TESTED_TYPE>(MATRIX_SIZE), get_omega(omega));
                en = std::chrono::high_resolution_clock::now();
                error += abs(answ - x);
                time += std::chrono::duration_cast<std::chrono::microseconds>(en-st).count();
            }
            file<<error/N_TO_AV<<","<<time/N_TO_AV<<",";             
        }
        error = 0;
        time = 0;
        for(int a=0; a < N_TO_AV; a++){
            st = std::chrono::high_resolution_clock::now();
            answ = SOR(A, b, NIter, 1e-7, std::vector<TESTED_TYPE>(MATRIX_SIZE));
            en = std::chrono::high_resolution_clock::now();
            error += abs(answ - x);
            time += std::chrono::duration_cast<std::chrono::microseconds>(en-st).count();
        }
        file<<error/N_TO_AV<<","<<time/N_TO_AV<<",";  
        file<<"\n";
    }
    file.close();
}