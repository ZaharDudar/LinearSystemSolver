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
    static Matrix<T> create_matrix_from_array(std::vector<T>, int, int);
    static Matrix<T> create_3diag_matrix(std::vector<T>,std::vector<T>,std::vector<T>);
};

//Create diagonal matrix from vector
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
        tmp_matrix(i+1,i,a[i]);
    }
    for(int i=0; i<tmp_matrix.nx-1; i++){
        tmp_matrix(i, i+1, c[i]);
    }
    return tmp_matrix;
}

// Create matrix from given array with given shape
template<typename T>
Matrix<T> Matrix<T>::create_matrix_from_array(std::vector<T> a, int nx, int ny){
    if(a.size() != (nx * ny)) {std::cout<<"Bad shape!"; throw 1;}

    Matrix<T> tmp_matrix;
    for(int i=0; i<nx*ny; i++){
        tmp_matrix.matrix = a;
    }
    tmp_matrix.nx = nx;
    tmp_matrix.ny = ny;
    return tmp_matrix;
}




// _________
template<typename T>
std::vector<T> SolveThomasAglorithm(Matrix<T> A, std::vector<T> B){
    std::vector<T> p;
    std::vector<T> q;
    std::vector<T> x;
    if(A.shape().first != A.shape().second){std::cout<<"Matrix not square!"; throw 1;}
    int n = A.shape().first;

    p.resize(n);
    q.resize(n);
    x.resize(n);

    //Прямой ход
    p[0] = - A(0,1) / A(0,0);
    q[0] = B[0] / A(0,0);

    for(int i=1; i<n-1; i++){
        p[i] = - A(i, i+1) / (A(i,i-1) * p[i-1] + A(i,i));
        q[i] = (B[i] - A(i,i-1) * q[i-1]) / (A(i,i-1) * p[i-1] + A(i, i));
    }

    //Обратный ход
    x[n-1] = (B[n-1] - A(n-1,n-2) * q[n-2]) / (A(n-1,n-2) * p[n-2] + A(n-1, n-1));
    for(int i=n-2; i>=0; i--){
        x[i] = x[i+1] * p[i] + q[i];
    }

    return x;
}

#endif 