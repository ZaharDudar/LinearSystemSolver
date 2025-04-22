#include <LinearSystemSolver.h>
#include <random>
#include <chrono>
#include <fstream>

#define TESTED_TYPE double

int main(){
    const int MATRIX_SIZE = 256;
    std::map<std::pair<int, int>, TESTED_TYPE> a;
    std::ifstream inputFile;
    inputFile.open("./256x256matrixDOK.txt");
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
    const int N_TO_AV=10;
    
    std::cout<<"load matrix\n";
    CSRMatrix<TESTED_TYPE> A(a,MATRIX_SIZE,MATRIX_SIZE);
    // CSRMatrix<TESTED_TYPE> A = CSRMatrix<TESTED_TYPE>::CSR_from_reg_matrix(Matrix<TESTED_TYPE>::create_matrix_from_array(a, MATRIX_SIZE, MATRIX_SIZE)); 
    std::cout<< reg_matrix_from_CSR(A);
    std::cout<< "max lambda= "<<FindMaxLambda(A,1e-14,500)<<"\n";
    std::cout<<"compute b\n";
    auto b = A*x;
    

    std::ofstream file;
    file.open("./UniversalBoost.csv");
    file<<"N_iter,Jacobi_error,Jacobi_time, Jacobi_Boosted_error, Jacobi_Boosted_time,G-S_error, G-S_time, G-S_Boosted_error, G-S_Boosted_time, Sym_G-S_error, Sym_G-S_time, Sym_G-S_Boosted_error, Sym_G-S_Boosted_time\n";
    auto st = std::chrono::high_resolution_clock::now();
    auto en = std::chrono::high_resolution_clock::now();
    double time;
    std::vector<TESTED_TYPE> answ;
    double error;
    std::cout<<"start writing\n";
    for(int NIter = 1; NIter <= 500; NIter+=1){
        if(NIter%10==0){std::cout<<((float)NIter/(float)500 * (float)100)<<"%\n";}
        // std::cout<<((float)NIter/(float)1000 * (float)100)<<"%\n";
        file<<NIter<<",";

        //Jacobi
        error = 0;
        time = 0;
        for(int a=0; a < N_TO_AV; a++){
            st = std::chrono::high_resolution_clock::now();
            answ = JacobiMethod(A, b, NIter, 1e-7, std::vector<TESTED_TYPE>(MATRIX_SIZE));
            en = std::chrono::high_resolution_clock::now();
            // std::cout<<answ<<"\n";
            error += abs(answ - x);
            time += std::chrono::duration_cast<std::chrono::microseconds>(en-st).count();
        }
        file<<error/N_TO_AV<<","<<time/N_TO_AV<<",";      


        error = 0;
        time = 0;
        for(int a=0; a < N_TO_AV; a++){
            auto step = [](CSRMatrix<TESTED_TYPE> &A, std::vector<TESTED_TYPE> &b, std::vector<TESTED_TYPE> &x0) { return JacobiMethod(A, b, 1, (TESTED_TYPE)1e-7, x0);};
            
            st = std::chrono::high_resolution_clock::now();
            answ = universalChebBoost<TESTED_TYPE>(step, A, b, NIter, (TESTED_TYPE)1e-7, std::vector<TESTED_TYPE>(MATRIX_SIZE));    
            en = std::chrono::high_resolution_clock::now();
            
            error += abs(answ - x);
            time += std::chrono::duration_cast<std::chrono::microseconds>(en-st).count();
        }
        file<<error/N_TO_AV<<","<<time/N_TO_AV<<",";

        //G-S
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
            auto step = [](CSRMatrix<TESTED_TYPE> &A, std::vector<TESTED_TYPE> &b, std::vector<TESTED_TYPE> &x0) { return GaussSeidelMethod(A, b, 1, (TESTED_TYPE)1e-7, x0);};
            
            st = std::chrono::high_resolution_clock::now();
            answ = universalChebBoost<TESTED_TYPE>(step, A, b, NIter, (TESTED_TYPE)1e-7, std::vector<TESTED_TYPE>(MATRIX_SIZE));    
            en = std::chrono::high_resolution_clock::now();
            
            error += abs(answ - x);
            time += std::chrono::duration_cast<std::chrono::microseconds>(en-st).count();
        }
        file<<error/N_TO_AV<<","<<time/N_TO_AV<<",";

        //Sym G-S
        error = 0;
        time = 0;
        for(int a=0; a < N_TO_AV; a++){
            st = std::chrono::high_resolution_clock::now();
            answ = SymmetricGSMethod(A, b, NIter, 1e-7, std::vector<TESTED_TYPE>(MATRIX_SIZE));
            en = std::chrono::high_resolution_clock::now();
            error += abs(answ - x);
            time += std::chrono::duration_cast<std::chrono::microseconds>(en-st).count();
        }
        file<<error/N_TO_AV<<","<<time/N_TO_AV<<",";      


        error = 0;
        time = 0;
        for(int a=0; a < N_TO_AV; a++){
            auto step = [](CSRMatrix<TESTED_TYPE> &A, std::vector<TESTED_TYPE> &b, std::vector<TESTED_TYPE> &x0) { return SymmetricGSMethod(A, b, 1, (TESTED_TYPE)1e-7, x0);};
            
            st = std::chrono::high_resolution_clock::now();
            answ = universalChebBoost<TESTED_TYPE>(step, A, b, NIter, (TESTED_TYPE)1e-7, std::vector<TESTED_TYPE>(MATRIX_SIZE));    
            en = std::chrono::high_resolution_clock::now();
            
            error += abs(answ - x);
            time += std::chrono::duration_cast<std::chrono::microseconds>(en-st).count();
        }
        file<<error/N_TO_AV<<","<<time/N_TO_AV<<"\n";

    }
    file.close();

}