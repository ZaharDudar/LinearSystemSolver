#pragma once

#include <vector>
#include <span>
#include <fstream>

template<typename T>
class ColMajMatrix
{
private:
    int nc;
    int nr;
    std::vector<T> matrix;
public:
    ColMajMatrix(){nc = 0; nr = 0;};
    ColMajMatrix(int nr, int nc){this->nc= nc; this->nr = nr; matrix.resize(nc*nr);};
    std::pair<int, int> shape() const { return std::make_pair(nr, nc);}
    // T& operator()(int i, int j){ return matrix[i * nc + j];} //Проблема если где-то прописать vector.push_back(mtrx(i,j)), в вектор кажется добавится ссылка на элемент и его можно будет поменять оттуда 
    T operator()(int i, int j) const { return matrix[i + j * nr];} 
    void operator()(int i, int j, T set_val){ matrix[i + j * nr]=set_val;}
    std::vector<T> operator*(const std::vector<T>&);
    ColMajMatrix<T> operator*(const ColMajMatrix<T>&);


    std::span<T> get_col(int j_col);
    std::span<T> get_col(int j_col, int st_id, int end_id);

    ColMajMatrix<T> get_transposed() const {
        ColMajMatrix<T> tmp(nc,nr);
        for(int j=0; j<nc; j++){
            for(int i=0; i<nr; i++){
                tmp(j, i, (*this)(i,j));
            }
        }
        return tmp;
    };

    void addCol(const std::vector<T>&);

    static ColMajMatrix<T> create_eye(int n);
    static ColMajMatrix<T> create_diag_matrix(const std::vector<T>&);
    static ColMajMatrix<T> create_matrix_from_array(const std::vector<T>&, int, int);
    static ColMajMatrix<T> create_3diag_matrix(const std::vector<T>&,const std::vector<T>&,const std::vector<T>&);
};

//Create diagonal matrix from vector
template<typename T>
ColMajMatrix<T> ColMajMatrix<T>::create_diag_matrix(const std::vector<T>& diag){
    ColMajMatrix<T> tmp_matrix;
    tmp_matrix.nc = diag.size();
    tmp_matrix.nr = diag.size();
    for(int i=0; i<tmp_matrix.nc * tmp_matrix.nr; i++){
        tmp_matrix.matrix.push_back(0);
    }
    for(int i=0; i<diag.size(); i++){
        tmp_matrix(i,i,diag[i]);
    }
    return tmp_matrix;
}

template<typename T>
ColMajMatrix<T> ColMajMatrix<T>::create_eye(int n){
    ColMajMatrix<T> tmp_matrix;
    tmp_matrix.nc = n;
    tmp_matrix.nr = n;
    for(int i=0; i<tmp_matrix.nc * tmp_matrix.nr; i++){
        tmp_matrix.matrix.push_back(0);
    }
    for(int i=0; i<n; i++){
        tmp_matrix(i,i,(T)1);
    }
    return tmp_matrix;
}



//Create 3-diagonal matrix with main diagonal b. Essential: lenght b = lenght a + 1 = lenght c + 1 
template<typename T>
ColMajMatrix<T> ColMajMatrix<T>::create_3diag_matrix(const std::vector<T>& a ,const std::vector<T>& b,const std::vector<T>& c){
    if(a.size() != c.size() or a.size() != (b.size()-1)) {std::cout<<"Bad diagonals lenghts!"; throw 1;}

    ColMajMatrix<T> tmp_matrix = ColMajMatrix<T>::create_diag_matrix(b);
    for(int i=0; i<tmp_matrix.nc-1; i++){
        tmp_matrix(i+1,i,a[i]);
    }
    for(int i=0; i<tmp_matrix.nc-1; i++){
        tmp_matrix(i, i+1, c[i]);
    }
    return tmp_matrix;
}

// Create matrix from given array with given shape
template<typename T>
ColMajMatrix<T> ColMajMatrix<T>::create_matrix_from_array(const std::vector<T>& a, int nr, int nc){
    if(a.size() != (nc * nr)) {std::cout<<"Bad shape!" <<a.size(); throw std::length_error(std::to_string(a.size()) + " != " + std::to_string(nc * nr));}

    ColMajMatrix<T> tmp_matrix;
    tmp_matrix.matrix.resize(nc*nr);
    for(int j=0; j<nc; j++){
        for(int i=0; i<nr; i++){
            tmp_matrix.matrix[j*nr+i] = a[j+i*nc];
        }
    }
    tmp_matrix.nc = nc;
    tmp_matrix.nr = nr;
    return tmp_matrix;
}

template<typename T>
std::vector<T> ColMajMatrix<T>::operator*(const std::vector<T>& B){
    if(this->shape().second != B.size()) {std::cout<<"lenght row != length B"; throw 1;}
    std::vector<T> result(this->shape().first);

    T tmp_sum=0;
    for(int j=0; j<this->shape().second; j++){
        tmp_sum = 0;
        for(int i=0; i<this->shape().first; i++){
            result[i] += (*this)(i,j) * B[j];
        }
    }
    return result;
}


template<typename T>
std::ostream& operator<<(std::ostream& str,const ColMajMatrix<T>& A){
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
std::span<T> ColMajMatrix<T>::get_col(int j_col){
    std::span<T> col(matrix.begin()+j_col * nr, nr);
    return col;
}
template<typename T>
std::span<T> ColMajMatrix<T>::get_col(int j_col, int st_id, int end_id){
    std::span<T> col(matrix.begin() + j_col * nr + st_id,end_id-st_id);
    return col;
}


template<typename T>
void ColMajMatrix<T>::addCol(const std::vector<T>& col){
    this->nc++;
    for(int i=0; i<col.size(); i++){
        this->matrix.push_back(col[i]);
    }

}