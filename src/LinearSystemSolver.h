#pragma once
#include <vector>
#include <iostream>


#include <Matrix.h>
#include <CSRMatrix.h>
#include <VectorMath.h>
#include <IterMethods.h>

//FIRST IND - ROW, SECOND IND - COLUMN


// _________
template<typename T>
std::vector<T> SolveThomasAglorithm(const Matrix<T>& A,const std::vector<T>& B){
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

template<typename T>
std::vector<T> SolveThomasAglorithm(const std::vector<T>& a, const std::vector<T>& b, const std::vector<T>& c, const std::vector<T>& B){
    std::vector<T> p;
    std::vector<T> q;
    std::vector<T> x;
    if((a.size() != c.size()) and b.size()-1 != a.size()){std::cout<<"Wrong lenghts!"; throw 1;}
    int n = b.size();

    p.resize(n);
    q.resize(n);
    x.resize(n);

    //Прямой ход
    p[0] = - c[0] / b[0];
    q[0] = B[0] / b[0];

    for(int i=1; i<n-1; i++){
        p[i] = - c[i] / (a[i-1] * p[i-1] + b[i]);
        q[i] = (B[i] - a[i-1] * q[i-1]) / (a[i-1] * p[i-1] + b[i]);
    }

    //Обратный ход
    x[n-1] = (B[n-1] - a[n-2] * q[n-2]) / (a[n-2] * p[n-2] + b[n-1]);
    for(int i=n-2; i>=0; i--){
        x[i] = x[i+1] * p[i] + q[i];
    }

    return x;
}
template<typename T>
std::pair<Matrix<T>, Matrix<T>> hausholder_alorithm(const Matrix<T>& A){
    //QR разлжения алгоритмом хаусхолдера

    int m = A.shape().first;
    int n = A.shape().second;

    Matrix<T> R = A;
    Matrix<T> Q = Matrix<T>::create_eye(m);

    for(int j=0; j < n; j++){
        std::vector<T> tmp = R.get_col(j, j, m);
        std::vector<T> e1 = tmp*(T)0;
        e1[0]=1;
        std::vector<T> v = tmp - abs(tmp)*e1;
        if(abs(v) == 0){continue;}

        for(int col = 0; col < n; col++){
            std::vector<T> x = R.get_col(col, j, m); 
            x = x - 2*(v*x)/(v*v)*v;
            for(int i = j; i < m; i++){
                R(i, col, x[i-j]);
            }
        }
        
        for(int row = 0; row < m; row++){
            std::vector<T> e = Q.get_row(row, j, m); 
            e = e - 2*(v*e)/(v*v)*v;
            
            for(int i = j; i < m; i++){
                Q(row, i, e[i-j]);
            }
        }
    }
    return std::make_pair(Q,R);
}

template<typename T>
std::vector<T> solve_by_QR(const Matrix<T>& Q, const Matrix<T>& R, const std::vector<T>& f){
    std::vector<T> b = Q.get_transposed() * f;
    Matrix<T> res = R;
    std::vector<T> x(f.size());
    

    for(int j=R.shape().second-1; j>=1; j--){
        for(int i = j-1; i>=0; i--){
            b[i] -= res(i,j) * b[j] / res(j,j); 
            // res(i,j, res(i,j) - res(i,j) * res(j,j) / res(j,j)); <- необязательные действия
        }
    }

    for(int i=0; i<f.size(); i++){
        x[i] = b[i]/res(i,i);
    }
    return x;
}
