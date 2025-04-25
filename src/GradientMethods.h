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

template<typename T>
std::vector<T> CG(const CSRMatrix<T>& A, const std::vector<T>& b, unsigned int nIter, T epsylon, const std::vector<T> x0){
    std::vector<T> x = x0;
    std::vector<T> r = A*x0 - b;
    std::vector<T> d = r;
    std::vector<T> r1 = r;

    for(unsigned int i=0; i<nIter; i++){
        x = x - r*r/(d*(A*d)) * d;
        r1 = A*x - b;
        d = r1 + r1*r1/(r*r)*d;
        r = r1;
        if(abs(r) <= epsylon){ return x;}
    }
    return x;
}

template<typename T>
std::vector<T> BiCG(const CSRMatrix<T>& A, const std::vector<T>& b, unsigned int nIter, T epsylon, const std::vector<T> x0){
    std::vector<T> x = x0;
    // std::vector<T> r1 = (T)(-1)*(A*x0-b);
    std::vector<T> r1 = A*x0-b;
    std::vector<T> d1 = r1;
    std::vector<T> r2 = r1;
    std::vector<T> d2 = r1;
    
    std::vector<T> rNew1 = r1;
    std::vector<T> rNew2 = r2;
    
    CSRMatrix<T> At = A.transpose();

    T alpha, betta;
    for(unsigned int i=0; i<nIter; i++){
        alpha = r2*r1/(d2 * (A*d1));
        x = x - alpha * d1;

        rNew1 = r1 - alpha * (A * d1);
        rNew2 = r2 - alpha * (At * d2);

        betta = rNew2 * rNew1 / (r2*r1);
        d1 = rNew1 + betta * d1;
        d2 = rNew1 + betta * d2;

        r1 = rNew1;
        r2 = rNew2;

        if(abs(r1) <= epsylon){ return x;}
    }
    return x;
}


template<typename T>
std::vector<T> CGS(const CSRMatrix<T>& A, const std::vector<T>& b, unsigned int nIter, T epsylon, const std::vector<T> x0){
    std::vector<T> x = x0;
    std::vector<T> R = A*x0-b;
    std::vector<T> r = R;
    std::vector<T> rNew = R;
    std::vector<T> d=R;
    std::vector<T> u=R;
    std::vector<T> q=R;

    T alpha, betta;
    for(unsigned int i=0; i<nIter; i++){
        alpha = R*r/(R * (A * d));
        q = u - alpha * (A * d);
        x = x - alpha * (u + q);
        
        rNew = r - alpha * (A * (u+q));

        betta = R*rNew/(R*r);
        u = rNew + betta*q;
        d = u + betta*(q + betta*d);

        r = rNew;
        if(abs(r) <= epsylon){ return x;}
    }
    return x;
}