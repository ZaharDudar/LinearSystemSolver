#pragma once

#include <vector>

template<typename T>
std::vector<T> operator+(std::vector<T> A, std::vector<T> B){
    if(A.size()!=B.size()) {std::cout<<"Different lenght of vectors!"; throw 1;}
    int n = A.size();

    std::vector<T> tmp(n);
    for(int i=0; i<n; i++){
        tmp[i] = A[i] + B[i];
    }
    return tmp;
}


template<typename T>
std::vector<T> operator-(std::vector<T> A, std::vector<T> B){
    if(A.size()!=B.size()) {std::cout<<"Different lenght of vectors!"; throw 1;}
    int n = A.size();

    std::vector<T> tmp(n);
    for(int i=0; i<n; i++){
        tmp[i] = A[i] - B[i];
    }
    return tmp;
}

template<typename T>
std::vector<T> operator*(std::vector<T> A, int B){
    int n = A.size();
    std::vector<T> tmp(n);
    for(int i=0; i<n; i++){
        tmp[i] = A[i] * B;
    }
    return tmp;
}

template<typename T>
std::vector<T> operator*(int A, std::vector<T> B){
    int n = B.size();
    std::vector<T> tmp(n);
    for(int i=0; i<n; i++){
        tmp[i] = A * B[i];
    }
    return tmp;
}

template<typename T>
T operator*(std::vector<T> A, std::vector<T> B){
    if(A.size()!=B.size()) {std::cout<<"Different lenght of vectors!"; throw 1;}
    int n = A.size();

    T tmp=0;
    for(int i=0; i<n; i++){
        tmp += A[i] * B[i];
    }
    return tmp;
}

