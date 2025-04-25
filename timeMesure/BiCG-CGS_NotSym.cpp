#include <LinearSystemSolver.h>
#include <random>
#include <chrono>
#include <fstream>

#define TESTED_TYPE double

int main(){
    const int MATRIX_SIZE = 256;
    std::map<std::pair<int, int>, TESTED_TYPE> a;
    std::ifstream inputFile;
    inputFile.open("./256x256matrixDOKNotSym.txt");
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
    file.open("./BiCG-CGS-NotSym.csv");
    file<<"N_iter,BiCG_error,BiCG_time,CSG_error,CSG_time\n";
    auto st = std::chrono::high_resolution_clock::now();
    auto en = std::chrono::high_resolution_clock::now();
    double time;
    std::vector<TESTED_TYPE> answ;
    double error;
    std::cout<<"start writing\n";
    for(int NIter = 1; NIter <= 500; NIter+=1){
        if(NIter%10==0){std::cout<<((float)NIter/(float)500 * (float)100)<<"%\n";}
        file<<NIter<<",";

        //BiCG
        error = 0;
        time = 0;
        for(int a=0; a < N_TO_AV; a++){
            st = std::chrono::high_resolution_clock::now();
            answ = BiCG(A, b, NIter, 1e-7, std::vector<TESTED_TYPE>(MATRIX_SIZE));
            en = std::chrono::high_resolution_clock::now();
            error += abs(answ - x);
            time += std::chrono::duration_cast<std::chrono::microseconds>(en-st).count();
        }
        file<<error/N_TO_AV<<","<<time/N_TO_AV<<",";      

        //CGS
        error = 0;
        time = 0;
        for(int a=0; a < N_TO_AV; a++){
            st = std::chrono::high_resolution_clock::now();
            answ = CGS(A, b, NIter, 1e-7, std::vector<TESTED_TYPE>(MATRIX_SIZE));
            en = std::chrono::high_resolution_clock::now();
            error += abs(answ - x);
            time += std::chrono::duration_cast<std::chrono::microseconds>(en-st).count();
        }
        file<<error/N_TO_AV<<","<<time/N_TO_AV<<"\n";      

    }
    file.close();

}