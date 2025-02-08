#include <LinearSystemSolver.h>
#include <random>
#include <chrono>
#include <fstream>

#define N_TO_AVERAGE 100

int main(){
    std::ofstream file;
    file.open("./TimeOfMulByVector.csv");
    file<<"N_non_zero, regular, CSR\n";

    int MATRIX_SIZE =1000;
    Matrix<int> regular_sparce_matrix(MATRIX_SIZE,MATRIX_SIZE);
    
    for(int i=0; i<MATRIX_SIZE; i++){
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            regular_sparce_matrix(i,j,0);
        }        
    }

    std::vector<int> B(MATRIX_SIZE);
    for(int i=0; i<MATRIX_SIZE; i++){
        B[i] = i+1;
    }
    long long int total_time=0;
    auto st = std::chrono::high_resolution_clock::now();
    auto en = std::chrono::high_resolution_clock::now();
    std::vector<int> res1;
    std::vector<int> res2;


    for(int N_non_zero = MATRIX_SIZE*MATRIX_SIZE/10; N_non_zero<=MATRIX_SIZE*MATRIX_SIZE; N_non_zero+=MATRIX_SIZE*MATRIX_SIZE/10){
        for(int i=0; i<MATRIX_SIZE; i++){
            for(int j=0; j<MATRIX_SIZE; j++){
                regular_sparce_matrix(i,j, 0);
            }
        }

        std::cout<<N_non_zero<<" N_non_zero \n";
        // for(int i=0; i<N_non_zero; i++){
        //     int ix = std::rand()%(MATRIX_SIZE-1);
        //     int jx = std::rand()%%(MATRIX_SIZE-1);
        //     regular_sparce_matrix(ix,jx, std::rand()%100+1);
        // }
        for(int i=0; i<MATRIX_SIZE; i++){
            for(int j=0; j<N_non_zero/MATRIX_SIZE; j++){
                regular_sparce_matrix(i,j, std::rand()%100+50);
            }
        }
        
        auto CSR_sparce_matrix = CSRMatrix<int>::CSR_from_reg_matrix(regular_sparce_matrix);
        std::cout<<"int CSR m has "<< CSR_sparce_matrix.get_N()<<" non zero els\n";
        std::cout<<"reg mtrx\n";
        total_time=0;
        for(int i=0; i<N_TO_AVERAGE; i++){
            st = std::chrono::high_resolution_clock::now();
            res1 = regular_sparce_matrix * B;
            en = std::chrono::high_resolution_clock::now();
            total_time += std::chrono::duration_cast<std::chrono::nanoseconds>(en-st).count();
        }
        total_time/=N_TO_AVERAGE;
        file << N_non_zero<<","<< total_time<<",";

        std::cout<<"CSR mtrx\n";
        total_time=0;
        for(int i=0; i<N_TO_AVERAGE; i++){
            st = std::chrono::high_resolution_clock::now();
            res2 = CSR_sparce_matrix * B;
            en = std::chrono::high_resolution_clock::now();
            total_time += std::chrono::duration_cast<std::chrono::nanoseconds>(en-st).count();
        }

        for(int i=0; i<res1.size(); i++){
            if(res1[i]!=res2[i]){
                throw 1;
            }
        }
        total_time/=N_TO_AVERAGE;
        file << total_time <<"\n";

    }

    file.close();
    return 0;
}