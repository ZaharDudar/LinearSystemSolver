#pragma once
#include <vector>


template<typename T>
class Matrix
{
private:
    int nc;
    int nr;
    std::vector<T> matrix;
public:
    Matrix(){nc = 0; nr = 0;};
    std::pair<int, int> shape(){ return std::make_pair(nr, nc);}
    // T& operator()(int i, int j){ return matrix[i * nc + j];} //Проблема если где-то прописать vector.push_back(mtrx(i,j)), в вектор кажется добавится ссылка на элемент и его можно будет поменять оттуда 
    T operator()(int i, int j){ return matrix[i * nc + j];} 
    T operator()(int i, int j, int set_val){ return matrix[i * nc + j]=set_val;}
    std::vector<T> operator*(std::vector<T>B);
    ~Matrix(){};
    static Matrix<T> create_diag_matrix(std::vector<T>);
    static Matrix<T> create_matrix_from_array(std::vector<T>, int, int);
    static Matrix<T> create_3diag_matrix(std::vector<T>,std::vector<T>,std::vector<T>);

};

//Create diagonal matrix from vector
template<typename T>
Matrix<T> Matrix<T>::create_diag_matrix(std::vector<T> diag){
    Matrix<T> tmp_matrix;
    tmp_matrix.nc = diag.size();
    tmp_matrix.nr = diag.size();
    for(int i=0; i<tmp_matrix.nc * tmp_matrix.nr; i++){
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
    for(int i=0; i<tmp_matrix.nc-1; i++){
        tmp_matrix(i+1,i,a[i]);
    }
    for(int i=0; i<tmp_matrix.nc-1; i++){
        tmp_matrix(i, i+1, c[i]);
    }
    return tmp_matrix;
}

// Create matrix from given array with given shape
template<typename T>
Matrix<T> Matrix<T>::create_matrix_from_array(std::vector<T> a, int nr, int nc){
    if(a.size() != (nc * nr)) {std::cout<<"Bad shape!"; throw 1;}

    Matrix<T> tmp_matrix;
    for(int i=0; i<nc*nr; i++){
        tmp_matrix.matrix = a;
    }
    tmp_matrix.nc = nc;
    tmp_matrix.nr = nr;
    return tmp_matrix;
}

template<typename T>
std::vector<T> Matrix<T>::operator*(std::vector<T>B){
    if(this->shape().second != B.size()) {std::cout<<"lenght row != length B"; throw 1;}
    std::vector<T> result(this->shape().first);

    T tmp_sum=0;
    for(int i=0; i<this->shape().first; i++){
        tmp_sum = 0;
        for(int j=0; j<this->shape().second; j++){
            tmp_sum += (*this)(i,j) * B[j];
        }
        result[i] = tmp_sum;
    }
    return result;
}