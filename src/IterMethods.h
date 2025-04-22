#pragma once

#include <Matrix.h>
#include <VectorMath.h>
#include <CSRMatrix.h>
// Here:
// Prime Iter Method, Jacobi Method, Gauss-Seidel Method, Symmetric GS Method, Chebyshev boosting PIM

template<typename T>
std::vector<T> PrimeIterMethod(const CSRMatrix<T>& A, const std::vector<T>& b, unsigned int nIter, T epsylon, const std::vector<T> x0){
    std::vector<T> X = x0;
    const T tau = 0.005;

    for(unsigned int i=0; i<nIter; i++){
        X = X - tau * (A * X - b); 
        if(abs(A * X - b) <= epsylon){ return X;}
    }
    return X;
}

template<typename T>
std::vector<T> JacobiMethod(const CSRMatrix<T>& A, const std::vector<T>& b, unsigned int nIter, T epsylon, const std::vector<T> x0){
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
        if(abs(A * X - b) <= epsylon){ return X;}
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
        if(abs(A * X - b) <= epsylon){ return X;}
    }
    return X;
}

template<typename T>
std::vector<T> SymmetricGSMethod(const CSRMatrix<T>& A, const std::vector<T>& b, unsigned int nIter, T epsylon, const std::vector<T> x0){
    std::vector<T> X = x0;
    std::vector<T> bUX = x0;
    std::vector<T> bLX = x0;
    int n = x0.size();

    T tmpSum=0;
    for(unsigned int iter=0; iter < nIter; iter++){
        //Iteration 1/2
        for(int i = 0; i < n; i++){
            tmpSum = 0;
            for(int j = i+1; j < n; j++){
                tmpSum += A(i, j) * X[j];
            }
            bUX[i] = b[i] - tmpSum;
        }
        X[0] = bUX[0]/A(0, 0);
        for(int i=1; i < n; i++){
            tmpSum = bUX[i];    
            for(int j=0; j<i; j++){
                tmpSum -= A(i,j) * X[j];
            }
            X[i] = tmpSum/A(i,i);
        }

        // Iteration 1
        for(int i = 0; i < n; i++){
            tmpSum = 0;
            for(int j = 0; j < i; j++){
                tmpSum += A(i, j) * X[j];
            }
            bLX[i] = b[i] - tmpSum;
        }
        X[n-1] = bLX[n-1]/A(n-1, n-1);
        for(int i=n-2; i >= 0; i--){
            tmpSum = bLX[i];    
            for(int j=i+1; j < n; j++){
                tmpSum -= A(i,j) * X[j];
            }
            X[i] = tmpSum/A(i,i);
        }
        if(abs(A * X - b) <= epsylon){ return X;}
    }
    return X;
}


template<typename T>
T FindMaxLambda(const CSRMatrix<T>& A, T epsylon, unsigned int nIter){
    std::vector<T> r = std::vector<T>(A.shape().first);
    r[0]=1;
    r[1] =9;
    T lambda;
    T Lastlambda;
    for(int iter=0; iter<nIter; iter++){
        r = A*r;
        r = r * (1 / abs(r));
        lambda = r * (A * r) / (r * r);
        if(std::abs(lambda-Lastlambda) <= epsylon) {return lambda;}
        Lastlambda=lambda;
    }
    return lambda;
}

template<typename T>
std::vector<T> SOR(const CSRMatrix<T>& A, const std::vector<T>& b, unsigned int nIter, T epsylon, const std::vector<T> x0){
    std::vector<T> X = x0;
    int n = x0.size();

    T firstSum = 0;
    T secondSum = 0;
    
    Matrix<T> DA(n,n);
    
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            if(i==j){
                DA(i, j, 1 - A(i,j) / A(i,i));
            }
            else{
                DA(i, j, A(i,j) / A(i,i));
            }
        }
    }
    
    T mu = FindMaxLambda(CSRMatrix<T>::CSR_from_reg_matrix(DA), epsylon, nIter);
    T omega = 1 + pow((mu / (1 + (T)sqrt(1-mu*mu))), 2);

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
            
            X[k] = (1-omega)*X[k] + omega*(b[k] - firstSum - secondSum)/A(k,k);
        }
        if(abs(A * X - b) <= epsylon){ return X;}
    }
    return X;
}

template<typename T>
std::vector<T> SOR(const CSRMatrix<T>& A, const std::vector<T>& b, unsigned int nIter, T epsylon, const std::vector<T> x0, T omega){
    std::vector<T> X = x0;
    int n = x0.size();

    T firstSum = 0;
    T secondSum = 0;
    
    Matrix<T> DA(n,n);
    
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            if(i==j){
                DA(i, j, 1 - A(i,j) / A(i,i));
            }
            else{
                DA(i, j, A(i,j) / A(i,i));
            }
        }
    }
    
    // T mu = FindMaxLambda(CSRMatrix<T>::CSR_from_reg_matrix(DA), epsylon, nIter);
    // T omega = 1 + pow((mu / (1 + (T)sqrt(1-mu*mu))), 2);

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
            
            X[k] = (1-omega)*X[k] + omega*(b[k] - firstSum - secondSum)/A(k,k);
        }
        if(abs(A * X - b) <= epsylon){ return X;}
    }
    return X;
}

std::vector<int> getPerms(int r, std::vector<int> in);


template<typename T>
std::vector<T> ChebSpeededPIM(const CSRMatrix<T>& A, const std::vector<T>& b, unsigned int nIter, T epsylon, const std::vector<T> x0){
    
    
    const int N_ROOTS = 64;
    std::vector<T> roots(N_ROOTS);
    T cosPiN = (T)cos(M_PI/N_ROOTS);
    T sinPiN = (T)sin(M_PI/N_ROOTS);
    roots[0] = (T)cos(M_PI/(2 * N_ROOTS));
    for(int i=1; i<N_ROOTS; i++){
        roots[i] = roots[i-1] * cosPiN - sqrt(1-roots[i-1]*roots[i-1]) * sinPiN;
    }
    T lambaMin = (T)1.1;
    T lambaMax = FindMaxLambda(A, epsylon, 1000);

    for(int i=0; i<N_ROOTS; i++){
        roots[i] = (lambaMin + lambaMax)/2 + (lambaMax - lambaMin)/2 * roots[i];
    }
    
    std::vector<int> perm = getPerms(6, {0,1});
    std::vector<T> X = x0;
    
    for(unsigned int i=0; i<nIter; i++){
        for(int j=0; j < N_ROOTS; j++){
            X = X - 1/roots[perm[j]] * (A * X - b); 
        }
        if(abs(A * X - b) <= epsylon){ return X;}
    }
    return X;
}
