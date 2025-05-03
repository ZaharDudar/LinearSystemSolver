#pragma once

#include <Matrix.h>
#include <ColMajMatrix.h>
#include <VectorMath.h>
#include <CSRMatrix.h>

//square
template<typename T>
struct UpperTriangularMatrix
{
    private:
        std::vector<T> data;
        int n=0;
    public:
        std::pair<int, int> shape() const { return std::make_pair(n, n);}
        T operator()(int i, int j){
            if(i>j) return 0;
            return data[(i+1)*i/2 + j];
        }
        void operator()(int i, int j, T val){
            data[(i+1)*i/2 + j];
        }
        void addCol(const std::vector<T>& col){
            n++;
            for(int i=0; i<n; i++){
                data.push_back(col[i]);
            }
        }
        void addCol(const std::span<T> col){
            n++;
            for(int i=0; i<col.size(); i++){
                data.push_back(col[i]);
            }
        }
};

template<typename T>
std::ostream& operator<<(std::ostream& str,struct UpperTriangularMatrix<T>& A){
    int m = A.shape().first;
    int n = A.shape().second;
    
    for(int i=0; i<m; i++){
        str << "[ ";
        for(int j=0; j<n; j++){
            str<<A(i,j)<<" ";
        }
        str<<"]\n";
    }
    return str;
} 


template<typename T>
struct Rotation
{
    double cos;
    double sin; 
    int i;
    std::vector<T> rotate(const std::vector<T>& v){
        std::vector<T> newVec = v;
        T tmp = newVec[0];
        newVec[0] = cos * tmp - sin * newVec[i];
        newVec[i] = sin * tmp + cos * newVec[i];
        return newVec;
    }
};

template<typename T>
struct Arnoldi
{
    int iteration=0;
    struct UpperTriangularMatrix<T> R;
    std::vector<Rotation<T>> rs;
    ColMajMatrix<T> basis;
    std::vector<T> applyAllRotations(const std::vector<T>& v){
        std::vector<T> res=v;
        for(int i=0; i<rs.size(); i++){
            res = rs[i].rotate(res);
        }
        return res;
    }
    void calcNewVect(const CSRMatrix<T>& A){
        //Считаем столбец матрицы Хессинберга и новый базисный вектор
        std::vector<T> v = A*basis.get_col(iteration);
        std::vector<T> h(iteration+2);

        for(int i=0; i<iteration+1; i++){
            h[i] = v * basis.get_col(i);
            v = v - h[i] * basis.get_col(iteration);
        }
        h[iteration+1] = abs(v);
        v = v * (1/h[iteration+1]);
        basis.addCol(v);

        //Добавляем вектор в R
        struct Rotation<T> newRot;
        h = applyAllRotations(h);
        newRot.cos = -h[0]/sqrt(h[0]*h[0] + h[iteration+1]*h[iteration+1]);
        newRot.sin = h[iteration+1]/sqrt(h[0]*h[0] + h[iteration+1]*h[iteration+1]);
        newRot.i = h.size()-1;
        rs.push_back(newRot);
        h = newRot.rotate(h);
        R.addCol(h);
        std::cout<<"R \n"<<R<<"\n";
        std::cout<<"Basis \n"<<basis<<"\n";
        iteration++;
    }
};

template<typename T>
std::vector<T> GMRES(const CSRMatrix<T>& A, const std::vector<T>& b, unsigned int nIter, T epsylon, const std::vector<T> x0){
    // std::vector<T>
    struct Arnoldi<T> arn;
    std::vector<T> r0 = A*x0-b;
    std::cout<<"r0: \n"<<r0<<"\n";

    T rho = abs(r0);
    arn.basis = ColMajMatrix<T>(x0.size(),1);
    for(int i=0; i<r0.size(); i++){
        arn.basis(i,0, r0[i]/rho);
    }
    std::cout<<"basis: \n"<<arn.basis<<"\n";
    T gamma = 0;

    for(int i=1; i<nIter+1; i++){
        arn.calcNewVect(A);
        std::vector<T> e1(i);
        e1[0]=1;
        gamma = arn.applyAllRotations(e1)[i-1];
        std::cout<<"gamma = "<<abs(gamma)<<"\n";
    }

    return r0;
}