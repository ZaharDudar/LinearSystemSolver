#pragma once

#include <Matrix.h>
#include <VectorMath.h>
#include <CSRMatrix.h>


template<typename T>
std::vector<T> SteepestGradientDescent(const CSRMatrix<T>& A, const std::vector<T>& b, unsigned int nIter, T epsylon, const std::vector<T> x0){
    std::vector<T> X = x0;
    T alpha;
    std::vector<T> r = A * X - b;

    for(unsigned int i=0; i<nIter; i++){
        r = A * X - b;
        alpha = r*r / (r*(A*r));
        X = X - alpha * r; 
        if(abs(r) <= epsylon){ return X;}
    }
    return X;
}
