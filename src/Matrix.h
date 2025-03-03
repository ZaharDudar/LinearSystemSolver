#pragma once
#include <vector>


template<typename T>
class Matrix
{
private:
    int nc;
    int nr;
    std::vector<T> matrix;
public:
    Matrix(){nc = 0; nr = 0;};
    Matrix(int nr, int nc){this->nc= nc; this->nr = nr; matrix.resize(nc*nr);};
    std::pair<int, int> shape() const { return std::make_pair(nr, nc);}
    // T& operator()(int i, int j){ return matrix[i * nc + j];} //Проблема если где-то прописать vector.push_back(mtrx(i,j)), в вектор кажется добавится ссылка на элемент и его можно будет поменять оттуда 
    T operator()(int i, int j) const { return matrix[i * nc + j];} 
    void operator()(int i, int j, T set_val){ matrix[i * nc + j]=set_val;}
    std::vector<T> operator*(const std::vector<T>&);
    Matrix<T> operator*(const Matrix<T>&);

    std::vector<T> get_col(int j_col);
    std::vector<T> get_col(int j_col, int st_id, int end_id);
    std::vector<T> get_row(int i_row);
    std::vector<T> get_row(int i_row, int st_id, int end_id);

    Matrix<T> get_transposed() const {
        Matrix<T> tmp(nc,nr);
        for(int i=0; i<nr; i++){
            for(int j=0; j<nc; j++){
                tmp(j, i, (*this)(i,j));
            }
        }
        return tmp;
    };

    static Matrix<T> create_eye(int n);
    static Matrix<T> create_diag_matrix(const std::vector<T>&);
    static Matrix<T> create_matrix_from_array(const std::vector<T>&, int, int);
    static Matrix<T> create_3diag_matrix(const std::vector<T>&,const std::vector<T>&,const std::vector<T>&);
};

//Create diagonal matrix from vector
template<typename T>
Matrix<T> Matrix<T>::create_diag_matrix(const std::vector<T>& diag){
    Matrix<T> tmp_matrix;
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
Matrix<T> Matrix<T>::create_eye(int n){
    Matrix<T> tmp_matrix;
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
Matrix<T> Matrix<T>::create_3diag_matrix(const std::vector<T>& a ,const std::vector<T>& b,const std::vector<T>& c){
    if(a.size() != c.size() or a.size() != (b.size()-1)) {std::cout<<"Bad diagonals lenghts!"; throw 1;}

    Matrix<T> tmp_matrix = Matrix<T>::create_diag_matrix(b);
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
Matrix<T> Matrix<T>::create_matrix_from_array(const std::vector<T>& a, int nr, int nc){
    if(a.size() != (nc * nr)) {std::cout<<"Bad shape!" <<a.size(); throw std::length_error(std::to_string(a.size()) + " != " + std::to_string(nc * nr));}

    Matrix<T> tmp_matrix;
    for(int i=0; i<nc*nr; i++){
        tmp_matrix.matrix = a;
    }
    tmp_matrix.nc = nc;
    tmp_matrix.nr = nr;
    return tmp_matrix;
}

template<typename T>
std::vector<T> Matrix<T>::operator*(const std::vector<T>& B){
    if(this->shape().second != B.size()) {std::cout<<"lenght row != length B"; throw 1;}
    std::vector<T> result(this->shape().first);

    T tmp_sum=0;
    for(int i=0; i<this->shape().first; i++){
        tmp_sum = 0;
        for(int j=0; j<this->shape().second; j++){
            tmp_sum += (*this)(i,j) * B[j];
        }
        result[i] = tmp_sum;
    }
    return result;
}

template<typename T>
std::ostream& operator<<(std::ostream& str,const Matrix<T>& A){
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
Matrix<T> Matrix<T>::operator*(const Matrix<T>& B){
    Matrix<T> result(this->nr, B.nc);

    T tmp_sum=0;
    for(int i=0; i<this->nr; i++){
        for(int j=0; j < B.nc; j++){
            tmp_sum = 0;
            for(int k=0; k < this->nc ; k++){
                tmp_sum += (*this)(i,k) * B(k,j);
            }
            result(i,j,tmp_sum);
        }
    }
    return result;
}

template<typename T>
std::vector<T> Matrix<T>::get_col(int j_col){
    std::vector<T> col(this->nr);
    for(int i=0; i<nr; i++){
        col[i] = (*this)(i, j_col);
    }
    return col;
}
template<typename T>
std::vector<T> Matrix<T>::get_col(int j_col, int st_id, int end_id){
    std::vector<T> col(end_id-st_id);
    for(int i=st_id; i<end_id; i++){
        col[i-st_id] = (*this)(i, j_col);
    }
    return col;
}

template<typename T>
std::vector<T> Matrix<T>::get_row(int i_row){
    std::vector<T> row(this->nc);
    for(int j=0; j<nc; j++){
        row[j] = (*this)(i_row, j);
    }
    return row;
}
template<typename T>
std::vector<T> Matrix<T>::get_row(int i_row, int st_id, int end_id){
    std::vector<T> row(end_id-st_id);
    for(int j=st_id; j<end_id; j++){
        row[j-st_id] = (*this)(i_row, j);
    }
    return row;
}
