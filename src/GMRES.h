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
        void clear(){
            n=0;
            data.clear();
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
            v = v - h[i] * basis.get_col(i);
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
std::vector<T> GMRES(const CSRMatrix<T>& A, const std::vector<T>& b, unsigned int nIter, unsigned int m, T epsylon, const std::vector<T> x0){
    // std::vector<T>
    struct Arnoldi<T> arn;
    std::vector<T> y = x0;
    std::vector<T> x = x0;
    std::vector<T> r0 = A*y-b;
    std::vector<T> e1(m);
    e1[0]=1;
    // std::cout<<"r0: \n"<<r0<<"\n";

    // T rho = abs(r0);
    // arn.basis = ColMajMatrix<T>(x0.size(),1);
    // for(int i=0; i<r0.size(); i++){
    //     arn.basis(i,0, r0[i]/rho);
    // }
    // std::cout<<"basis: \n"<<arn.basis<<"\n";
    T gamma = 0;
    T rho = 0;

    for(int i=0; i<2; i++){
        //Init Arnoldi
        r0 = A*x-b;
        e1 = e1*(T)0;
        e1[0]=1;

        rho = abs(r0);
        arn.basis = ColMajMatrix<T>(x0.size(),1);
        arn.R.clear();
        arn.rs.clear();
        arn.iteration=0;
        for(int i=0; i<r0.size(); i++){
            arn.basis(i,0, r0[i]/rho);
        }
        std::cout<<"rho = "<<rho<<"\n";

        //m iters of Arnoldi
        for(int j=0; j<m; j++){
            arn.calcNewVect(A);
            e1 = arn.rs[arn.iteration-1].rotate(e1);
            gamma = abs(e1[arn.iteration]);
            std::cout<<"gamma = "<<gamma<<"\n";
            if(abs(gamma)<=epsylon){ break;}
        }       
        std::cout<<arn.basis.get_col(0) * arn.basis.get_col(1)<<" <-dot of v\n"  ;
        std::cout<<"e1 = \n"<<e1<<"\n";
        //Ry = z
        y.clear();
        y.resize(arn.iteration,(T)0);

        y[arn.iteration-1] = e1[arn.iteration-1] / arn.R(arn.iteration-1,arn.iteration-1);
        for(int n=arn.iteration-2; n>=0; n--){
            y[n] = e1[n];
            for(int k=n+1; k<=arn.iteration-1; k++){
                y[n] -= arn.R(n, k) * y[k];
            }
            y[n]/=arn.R(n,n);
        }
        std::cout<<"y = \n"<<y<<"\n";
        x = x0 - arn.basis.partly_mul(y); 
        std::cout<<"x = \n"<<x<<"\n";
        if(abs(A * x - b) <= epsylon){ return x;}
    }

    return x;
}