#ifndef LSS_H
#define LSS_H

#include <vector>
#include <iostream>


#include <Matrix.h>
#include <SCRMatrix.h>


//FIRST IND - ROW, SECOND IND - COLUMN


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