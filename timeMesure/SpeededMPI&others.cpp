#include <LinearSystemSolver.h>
#include <random>
#include <chrono>
#include <fstream>

#define TESTED_TYPE double

int main(){
    std::vector<TESTED_TYPE> aOld = {
            20,3,2,0,3,0,0,2,3,0,3,22,3,3,0,0,0,0,3,3,2,3,20,0,1,0,0,0,4,0,0,3,0,21,0,4,0,1,0,2,3,0,1,0,22,3,0,2,3,3,0,0,0,4,3,20,4,0,3,0,0,0,0,0,0,4,20,4,1,0,2,0,0,1,2,0,4,20,0,0,3,3,4,0,3,3,1,0,22,0,0,3,0,2,3,0,0,0,0,23
        
        };
        
    const int MATRIX_SIZE = 128;
    const int N_TO_AV=5;
    std::map<std::pair<int, int>, TESTED_TYPE> a;
    std::ifstream inputFile;
    inputFile.open("./128x128matrixDOK.txt");
    std::string str;
    std::string tmpStr;
    std::cout<<"read file\n";
    int iIn,jIn, tmpi;
    TESTED_TYPE valIn;
    while (std::getline(inputFile, str))
    {
        for(int i=0; i< str.length(); i++){
            if(str[i] == ','){
                iIn = std::stoi(tmpStr);
                tmpi=i;
                tmpStr.clear();
                break;
            }
            tmpStr+=str[i];
        }
        tmpStr.clear();
        for(int i=tmpi+1; i< str.length(); i++){
            if(str[i] == ','){
                jIn = std::stoi(tmpStr);
                tmpi=i;
                tmpStr.clear();
                break;
            }
            tmpStr+=str[i];
        }
        tmpStr.clear();
        
        for(int i=tmpi+1; i< str.length(); i++){
            tmpStr+=str[i];
        }
        valIn = std::stod(tmpStr);
        a[std::make_pair(iIn, jIn)] = valIn;
        tmpStr.clear();
    }  
    inputFile.close();

    std::vector<TESTED_TYPE> x(MATRIX_SIZE);
    for(int i=0;i<MATRIX_SIZE;i++){
        x[i]=i+1;
    }

    std::cout<<"load matrix\n";
    CSRMatrix<TESTED_TYPE> A(a,MATRIX_SIZE,MATRIX_SIZE);
    // CSRMatrix<TESTED_TYPE> A = CSRMatrix<TESTED_TYPE>::CSR_from_reg_matrix(Matrix<TESTED_TYPE>::create_matrix_from_array(aOld, MATRIX_SIZE, MATRIX_SIZE)); 
    std::cout<< reg_matrix_from_CSR(A);
    std::cout<<"compute b\n";
    // CSRMatrix<TESTED_TYPE> A = CSRMatrix<TESTED_TYPE>::CSR_from_reg_matrix(Matrix<TESTED_TYPE>::create_matrix_from_array(aOld, MATRIX_SIZE, MATRIX_SIZE)); 
    auto b = A*x;
    
    std::ofstream file;
    file.open("./TimeMethods.csv");
    file<<"N_iter, MPI_error, MPI_time, Speded_MPI_error, Speded_MPI_time, Yakobi_error, Yakobi_time, GaussSeidel_error, GaussSeidel_time";
    file<<",SOR_error, SOR_time, SteepestGradientDescent_error, SteepestGradientDescent_time, SymmetricGSMethod_error, SymmetricGSMethod_time\n";
    auto st = std::chrono::high_resolution_clock::now();
    auto en = std::chrono::high_resolution_clock::now();
    double time;
    std::vector<TESTED_TYPE> answ;
    double error;
    for(int NIter = 1; NIter <= 1000; NIter+=1){
        if(NIter%10==0){std::cout<<((float)NIter/(float)1000 * (float)100)<<"%\n";}
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
            // for(int i=0; i<answ.size(); i++) {if(answ[i] - x[i] >= 1e-6) return 1;}
            en = std::chrono::high_resolution_clock::now();
            error += abs(answ - x);
            time += std::chrono::duration_cast<std::chrono::microseconds>(en-st).count();
        }
        file<<error/N_TO_AV<<","<<time/N_TO_AV<<",";
        
        error = 0;
        time = 0;
        for(int a=0; a < N_TO_AV; a++){
            st = std::chrono::high_resolution_clock::now();
            answ = JacobiMethod(A, b, NIter, 1e-7, std::vector<TESTED_TYPE>(MATRIX_SIZE));
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
        file<<error/N_TO_AV<<","<<time/N_TO_AV<<",";
        
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


        error = 0;
        time = 0;
        for(int a=0; a < N_TO_AV; a++){
            st = std::chrono::high_resolution_clock::now();
            answ = SteepestGradientDescent(A, b, NIter, 1e-7, std::vector<TESTED_TYPE>(MATRIX_SIZE));
            en = std::chrono::high_resolution_clock::now();
            error += abs(answ - x);
            time += std::chrono::duration_cast<std::chrono::microseconds>(en-st).count();
        }
        file<<error/N_TO_AV<<","<<time/N_TO_AV<<",";      


        error = 0;
        time = 0;
        for(int a=0; a < N_TO_AV; a++){
            st = std::chrono::high_resolution_clock::now();
            answ = SymmetricGSMethod(A, b, NIter, 1e-7, std::vector<TESTED_TYPE>(MATRIX_SIZE));
            en = std::chrono::high_resolution_clock::now();
            error += abs(answ - x);
            time += std::chrono::duration_cast<std::chrono::microseconds>(en-st).count();
        }
        file<<error/N_TO_AV<<","<<time/N_TO_AV<<"\n";        
    }
    file.close();

}