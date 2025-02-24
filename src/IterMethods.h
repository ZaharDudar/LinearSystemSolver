#pragma once

#include <Matrix.h>
#include <VectorMath.h>
#include <CSRMatrix.h>


template<typename T>
std::vector<T> PrimeIterMethod(const CSRMatrix<T>& A, const std::vector<T>& b, unsigned int nIter, T epsylon, const std::vector<T> x0){
    std::vector<T> X = x0;
    const T tau = 0.005;

    for(unsigned int i=0; i<nIter; i++){
        X = X - tau * (A * X - b); 
        if(abs(A * X - b) < epsylon){ return X;}
    }
    return X;
}

template<typename T>
std::vector<T> YakobyMethod(const CSRMatrix<T>& A, const std::vector<T>& b, unsigned int nIter, T epsylon, const std::vector<T> x0){
    std::vector<T> X = x0;

    std::vector<T> val;
    std::vector<int> cols;
    std::vector<int> rows;
    std::tie(val, cols, rows) = A.getRawData();

    T diagEl;
    T tmp_sum=0;
    std::vector<T> result(A.shape().first);
    for(unsigned int iter=0; iter < nIter; iter++){
        for(int i=0; i<A.shape().first; i++){
            tmp_sum = 0;
            diagEl = 0;
            for(int j=rows[i]; j<rows[i+1]; j++){
                if(i==cols[j]){
                    diagEl = val[j];
                }
                else{
                    tmp_sum += val[j] * X[cols[j]];
                }
            }
            result[i] = (b[i] - tmp_sum)/diagEl;
        }
        X = result;
        if(abs(A * X - b) < epsylon){ return X;}
    }
    return X;
}

template<typename T>
std::vector<T> GaussSeidelMethod(const CSRMatrix<T>& A, const std::vector<T>& b, unsigned int nIter, T epsylon, const std::vector<T> x0){
    std::vector<T> X = x0;
    int n = x0.size();

    T firstSum=0;
    T secondSum=0;
    for(unsigned int iter=0; iter < nIter; iter++){
        for(int k=0; k<n; k++){
            firstSum = 0;
            secondSum = 0;

            for(int j = k+1; j < n; j++){
                firstSum += A(k, j) * X[j];
            }
            for(int j = 0; j < k; j++){
                secondSum += A(k, j) * X[j];
            }
            
            X[k] = (b[k] - firstSum - secondSum)/A(k,k);
        }
        if(abs(A * X - b) < epsylon){ return X;}
    }
    return X;
}