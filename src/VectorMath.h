#pragma once

#include <cmath>
#include <vector>

template<typename T>
std::vector<T> operator+(const std::vector<T>& A,const std::vector<T>& B){
    // if(A.size()!=B.size()) {std::cout<<"Different lenght of vectors!"; throw 1;}
    int n = A.size();

    std::vector<T> tmp(n);
    for(int i=0; i<n; i++){
        tmp[i] = A[i] + B[i];
    }
    return tmp;
}


template<typename T>
std::vector<T> operator-(const std::vector<T>& A,const std::vector<T>& B){
    // if(A.size()!=B.size()) {std::cout<<"Different lenght of vectors!"; throw 1;}
    int n = A.size();

    std::vector<T> tmp(n);
    for(int i=0; i<n; i++){
        tmp[i] = A[i] - B[i];
    }
    return tmp;
}

template<typename T>
std::vector<T> operator*(const std::vector<T>& A, T B){
    int n = A.size();
    std::vector<T> tmp(n);
    for(int i=0; i<n; i++){
        tmp[i] = A[i] * B;
    }
    return tmp;
}

template<typename T>
std::vector<T> operator*(T A, const std::vector<T>& B){
    int n = B.size();
    std::vector<T> tmp(n);
    for(int i=0; i<n; i++){
        tmp[i] = A * B[i];
    }
    return tmp;
}

template<typename T>
std::vector<T> operator*(T A, const std::span<T>& B){
    int n = B.size();
    std::vector<T> tmp(n);
    for(int i=0; i<n; i++){
        tmp[i] = A * B[i];
    }
    return tmp;
}

template<typename T>
T operator*(const std::vector<T>& A,const std::vector<T>& B){
    if(A.size()!=B.size()) {std::cout<<"Different lenght of vectors!"; throw 1;}
    int n = A.size();

    T tmp=0;
    for(int i=0; i<n; i++){
        tmp += A[i] * B[i];
    }
    return tmp;
}

template<typename T>
T operator*(const std::vector<T>& A,const std::span<T>& B){
    if(A.size()!=B.size()) {std::cout<<"Different lenght of vectors!"; throw 1;}
    int n = A.size();

    T tmp=0;
    for(int i=0; i<n; i++){
        tmp += A[i] * B[i];
    }
    return tmp;
}

template<typename T>
T operator*(const std::span<T>& A,const std::vector<T>& B){
    if(A.size()!=B.size()) {std::cout<<"Different lenght of vectors!"; throw 1;}
    int n = A.size();

    T tmp=0;
    for(int i=0; i<n; i++){
        tmp += A[i] * B[i];
    }
    return tmp;
}

template<typename T>
T operator*(const std::span<T>& A,const std::span<T>& B){
    if(A.size()!=B.size()) {std::cout<<"Different lenght of vectors!"; throw 1;}
    int n = A.size();

    T tmp=0;
    for(int i=0; i<n; i++){
        tmp += A[i] * B[i];
    }
    return tmp;
}

template<typename T>
T abs(std::vector<T> A){
    return std::sqrt(A*A);
}

template<typename T>
std::ostream& operator<<(std::ostream& str, const std::vector<T>& A){
    for(int i=0; i<A.size(); i++){
        str<<"[ "<<A[i]<<" ]\n";
    }
    return str;
}

