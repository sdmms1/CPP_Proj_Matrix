//
// Created by Lenovo on 2020/6/8.
//

#ifndef PROJECT_N_MATRIX_HPP
#define PROJECT_N_MATRIX_HPP

#include <iostream>
#include "S_Matrix.hpp"
using namespace std;

struct SizeException:public exception{
    int type;
    explicit SizeException(int t){type = t;}
    const char * what() const throw(){
        if (type == 1) return "The size of the two matrices don't support the add operation";
        else if (type == 2) return "The size of the two matrices don't support the subtract operation";
        else if (type == 3) return "The size of the two matrices don't support the element-wise multiply operation";
        else if (type == 5) return "The cross product in n-dimension space should provide n-1 vectors";
        else return "The size of the two matrices don't support the multiply operation";
    }
};


template <class T>
class N_Matrix {
public:
    int col{};
    int row{};
    T * matrix;

    N_Matrix(){
        this->col = 1;
        this->row = 1;
    }
    N_Matrix(int row,int col){
        this->row = row;
        this->col = col;
        matrix = new T[col*row];
    }
    N_Matrix(int row,int col,T*Mat){
        this->col = col;
        this->row = row;
        matrix = new T[row*col];
        T *p = matrix;
        for (int i = 0; i < row*col; ++i)
            *p++ = *Mat++;

    }
    N_Matrix(S_Matrix<T>& other){
        this->row = other.row;
        this->col = other.col;
        matrix = new T [row*col];
        for (int i = 0; i < row*col; ++i) {
            matrix[i] = 0;
        }

        for (int i = 0; i < row; ++i) {
            while (other.row_offsets[i] == other.row_offsets[i+1]) i++;
            if (i>row) break;
            for (int j = other.row_offsets[i]; j < other.row_offsets[i+1]; ++j) {
                matrix[i*col+other.col_idx[j]] = other.values[j];
            }
        }

    }
    void setMatrix(T*Mat){
        T *p = matrix;
        for (int i = 0; i < row*col; ++i)
            *p++ = *Mat++;
    }
    void modifyVal(int row_num,int col_num,T val){
        matrix[row_num*col+col_num] = val;
    }

    template <typename T1>
    N_Matrix operator + (const N_Matrix<T1>& other) const{
        try{
            if ((this->col != other.col) | (this->row != other.row)) throw SizeException(1);
            T * result = new T[col*row];
            for (int i = 0; i < col*row; ++i) {
                result[i] = matrix[i] + other.matrix[i];
            }
            return N_Matrix{row,col,result};
        }catch (SizeException& e){
            cout<<e.what();
        }
    }

    template <typename T1>
    N_Matrix operator - (const N_Matrix<T1>& other) const{
        try{
            if ((this->col != other.col) | (this->row != other.row)) throw SizeException(2);
            T * result = new T[col*row];
            for (int i = 0; i < col*row; ++i) {
                result[i] = matrix[i] - other.matrix[i];
            }
            return N_Matrix{row,col,result};
        }catch (SizeException& e){
            cout<<e.what();
        }
    }

    template <typename T1>
    N_Matrix scalar_product (T1 val) const{
        T* result = new T[col*row];
        for (int i = 0; i < col*row; ++i) {
            result[i] = val*matrix[i];
        }
        return N_Matrix {row,col,result};
    }

    template <typename T1>
    N_Matrix operator / (T1 val){
        T* result = new T[col*row];
        for (int i = 0; i < col*row; ++i) {
            result[i] = matrix[i]/val;
        }
        return N_Matrix {row,col,result};
    }

    N_Matrix transposition(){
        T * result = new T[col*row];
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < col; ++j) {
                result[j*row+i] = matrix[i*col+j];
            }
        }
        return N_Matrix{col,row,result};
    }

    N_Matrix conjugation(){
        T * result = new T[col*row];
        for (int i = 0; i < col*row; ++i) {
            result[i] = conj(matrix[i]);
        }
        return N_Matrix {row,col,result};
    }

    template <typename T1>
    N_Matrix element_wise_mult(const N_Matrix<T1>& matrix2){
        try{
            if ((row != matrix2.row)|(col != matrix2.col)) throw SizeException(3);
            T * result = new T [col*row];
            for (int i = 0; i < row*col; ++i) {
                result[i] = matrix[i] * matrix2.matrix[i];
            }
            return N_Matrix{row,col,result};
        }catch (SizeException& e){
            cout<<e.what();
        }
    }

    template <typename T1>
    N_Matrix operator * (const N_Matrix<T1>& other) const {
        try{
            if (col != other.row) throw SizeException(4);
            T * result = new T [other.col*row];
            for (int i = 0; i < row; ++i) {
                for (int j = 0; j < other.col; ++j) {
                    T res = 0;
                    for (int k = 0; k < col; ++k) {
                        res = res+ matrix[i*col+k] * other.matrix[j+k*other.col];
                    }
                    result[i*other.col+j] = res;
                }
            }
            return N_Matrix{row,other.col,result};
        }catch (SizeException& e){
            cout<<e.what();
        }
    }

    template <typename T1>
    N_Matrix dot_product(const N_Matrix& other){
        try {
            if ((col != other.col) | (row != other.row)) throw SizeException(4);
            N_Matrix trans = this->transposition();
            return trans * other;
        }catch (SizeException& e){
            cout<<e.what();
        }
    }

    N_Matrix rotate(){
        T * result = new T [row*col];
        for (int i = 0; i < row*col; ++i) {
            result[i] = matrix[col*row-1-i];
        }
        return N_Matrix{row,col,result};
    }

    template <typename T1>
    N_Matrix convolution(N_Matrix<T1> & core){
        int row_cnt = row+2*core.row-2;
        int col_cnt = col+2*core.col-2;
        T * augment = new T [row_cnt*col_cnt];
        //rotate 180
        N_Matrix rot = core.rotate();
        //make up 0
        for (int i = 0; i < row_cnt; ++i) {
            for (int j = 0; j < col_cnt; ++j) {
                if (i<core.row-1 | i>=core.row-1+row | j<core.col-1 | j>=core.col-1+col) augment[i * col_cnt + j] = 0;
                else augment[i * col_cnt + j] = matrix[(i - core.row + 1) * col + (j - core.col + 1)];
            }
        }
        T * result = new T [(row+core.row-1)*(col+core.col-1)];
        for (int i = 0; i < row+core.row-1; ++i) {
            for (int j = 0; j < col+core.col-1; ++j) {
                T res = 0;
                for (int k = 0; k < core.row; ++k) {
                    for (int l = 0; l < core.col; ++l) {
                        res += rot.matrix[k*core.col+l]*augment[(i+k)*col_cnt+j+l];
                    }
                }
                result[i*(row+core.row-1)+j] = res;
            }
        }
        delete augment;
        return N_Matrix{row+core.row-1,col+core.col-1,result};
    }

    void printMatrix(){
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < col; ++j) {
                cout<<matrix[i*col+j]<<" ";
            }
            cout<<endl;
        }
    }
};



#endif //PROJECT_N_MATRIX_HPP