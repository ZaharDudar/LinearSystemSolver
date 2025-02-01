#ifndef LSS_H
#define LSS_H

#include <vector>
#include <iostream>


template<typename T>
class Matrix
{
private:
    int nx;
    int ny;
    std::vector<T> matrix;
public:
    Matrix(){nx = 0; ny = 0;};
    std::pair<int, int> shape(){ return std::make_pair(nx, ny);}
    T operator()(int i, int j){ return matrix[i * nx + j];}
    T operator()(int i, int j, int set_val){ return matrix[i * nx + j]=set_val;}
    ~Matrix(){};
    static Matrix<T> create_diag_matrix(std::vector<T>);
    static Matrix<T> create_3diag_matrix(std::vector<T>,std::vector<T>,std::vector<T>);
};

template<typename T>
Matrix<T> Matrix<T>::create_diag_matrix(std::vector<T> diag){
    Matrix<T> tmp_matrix;
    tmp_matrix.nx = diag.size();
    tmp_matrix.ny = diag.size();
    for(int i=0; i<tmp_matrix.nx * tmp_matrix.ny; i++){
        tmp_matrix.matrix.push_back(0);
    }
    for(int i=0; i<diag.size(); i++){
        tmp_matrix(i,i,diag[i]);
    }
    return tmp_matrix;
}


//Create 3-diagonal matrix with main diagonal b. Essential: lenght b = lenght a + 1 = lenght c + 1 
template<typename T>
Matrix<T> Matrix<T>::create_3diag_matrix(std::vector<T> a , std::vector<T> b, std::vector<T> c){
    if(a.size() != c.size() or a.size() != (b.size()-1)) {std::cout<<"Bad diagonals lenghts!"; throw 1;}

    Matrix<T> tmp_matrix = Matrix<T>::create_diag_matrix(b);
    for(int i=0; i<tmp_matrix.nx-1; i++){
        tmp_matrix(i,i+1,a[i]);
    }
    for(int i=0; i<tmp_matrix.nx-1; i++){
        tmp_matrix(i+1, i, c[i]);
    }
    return tmp_matrix;
}

#endif 