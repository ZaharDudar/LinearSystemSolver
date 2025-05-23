#pragma once

#include <vector>
#include <map>
#include <algorithm>
#include <Matrix.h>
#include <iostream>


template<typename T>
class CSRMatrix
{
private:
    int nc;
    int nr;
    int N_non_zero;
    std::vector<T> val;
    std::vector<int> cols;
    std::vector<int> rows;
    CSRMatrix(){};

public:
    CSRMatrix(const std::map<std::pair<int, int>, T>&,int, int);
    int get_N(){return N_non_zero;}
    std::pair<int, int> shape() const { return std::make_pair(nr, nc);}
    T operator()(int r, int c) const {
        for(int ind=rows[r]; ind<rows[r+1]; ind++){
            if(cols[ind] == c){
                return val[ind];
            }
        }
        return (T)0;
    }
    std::vector<T> operator*(const std::vector<T>&) const;
    std::vector<T> operator*(const std::span<T>&) const;
    std::tuple<std::vector<T>, std::vector<int>, std::vector<int>> getRawData() const;
    static CSRMatrix<T> CSR_from_reg_matrix(const Matrix<T>&);

    //Slow & Ugly.
    CSRMatrix<T> transpose() const;
};




template<typename T>
CSRMatrix<T>::CSRMatrix(const std::map<std::pair<int, int>, T>& in, int nc, int nr){
    this->nc=nc;
    this->nr=nr;
    N_non_zero=0;
    rows.push_back(0);

    // for (auto const &el : in){ 
    //     if(el.first.first >= nc){
    //         nc = el.first.first;
    //     }
    //     if(el.first.second >= nr){
    //         nr = el.first.second;
    //     }
    // }
    // //nr and nc after cycle is maximum indxs but it has to be number in a col and in a row relative
    // nr++; 
    // nc++

    // std::vector<T> zero_diag(nc*nr);
    // for(int i=0; i<nc*nr; i++) {zero_diag[i]=0;}

    // Matrix<T> tmp_matrix = Matrix<T>::create_matrix_from_array(zero_diag, nc, nr);
    // for (auto const &el : in){
    //     if(el.first.first >= nc or el.first.second >= nr){std::cout<<"Wrong shape!"; throw 1;}
    //     tmp_matrix(el.first.first, el.first.second, el.second);
    // }
    
    
    for(int i=0; i<nr; i++){
        for(int j=0; j<nc; j++){
            if(in.count(std::make_pair(i,j))!=0){
                val.push_back(in.at(std::make_pair(i,j)));
                cols.push_back(j);
                N_non_zero++;
            }
        }
        rows.push_back(N_non_zero);
    }

}

template<typename T>
Matrix<T> reg_matrix_from_CSR(const CSRMatrix<T>& A){
    int nr = A.shape().first;
    int nc = A.shape().second;

    Matrix<T> tmp_matrix(nr, nc);

    for(int i=0; i<nr; i++){
        for(int j=0; j<nr; j++){
            tmp_matrix(i,j,A(i,j));
        }
    }
    return tmp_matrix;
}


template<typename T>
CSRMatrix<T> CSRMatrix<T>::CSR_from_reg_matrix(const Matrix<T>& A){
    int nr = A.shape().first;
    int nc = A.shape().second;

    CSRMatrix<T> tmp_matrix;
    tmp_matrix.nc = nc;
    tmp_matrix.nr = nr;
    tmp_matrix.N_non_zero = 0;
    tmp_matrix.rows.push_back(0);

    for(int i=0; i<nr; i++){
        for(int j=0; j<nc; j++){
            if(A(i,j)!=0){
                tmp_matrix.val.push_back(A(i,j));
                tmp_matrix.cols.push_back(j);
                tmp_matrix.N_non_zero++;
            }
        }
        tmp_matrix.rows.push_back(tmp_matrix.N_non_zero);
    }
    return tmp_matrix;
}

template<typename T>
std::vector<T> CSRMatrix<T>::operator*(const std::vector<T>& B) const{
    std::vector<T> result(this->shape().first);

    T tmp_sum=0;
    for(int i=0; i<this->shape().first; i++){
        tmp_sum = 0;
        for(int j=this->rows[i]; j<this->rows[i+1]; j++){
            tmp_sum += this->val[j] * B[cols[j]];
        }
        result[i] = tmp_sum;
    }
    return result;
}

template<typename T>
std::vector<T> CSRMatrix<T>::operator*(const std::span<T>& B) const{
    std::vector<T> result(this->shape().first);

    T tmp_sum=0;
    for(int i=0; i<this->shape().first; i++){
        tmp_sum = 0;
        for(int j=this->rows[i]; j<this->rows[i+1]; j++){
            tmp_sum += this->val[j] * B[cols[j]];
        }
        result[i] = tmp_sum;
    }
    return result;
}

template<typename T>
std::tuple<std::vector<T>, std::vector<int>, std::vector<int>> CSRMatrix<T>::getRawData() const {
    return std::make_tuple(val, cols, rows);
}

template<typename T>
CSRMatrix<T> CSRMatrix<T>::transpose() const {
    Matrix<T> tmp(this->nr, this->nc);
    for(unsigned int i=0; i<this->nr; i++){
        for(unsigned int j=0; j<this->nr; j++){
            tmp(j, i, (*this)(i,j));
        }
    }
    return CSR_from_reg_matrix(tmp);
}