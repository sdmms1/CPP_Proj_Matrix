//
// Created by Lenovo on 2020/6/9.
//

#ifndef PROJECT_S_MATRIX_HPP
#define PROJECT_S_MATRIX_HPP

#include <iostream>
#include "N_Matrix.hpp"

using namespace std;


template<typename T>
class S_Matrix {
private:
    int col{};
    int row{};
    T *values;
    int *col_idx;
    int *row_offsets;
public:
    int getCol() const {
        return col;
    }

    void setCol(int col) {
        S_Matrix::col = col;
    }

    int getRow() const {
        return row;
    }

    void setRow(int row) {
        S_Matrix::row = row;
    }

    T *getValues() const {
        return values;
    }

    void setValues(T *values) {
        S_Matrix::values = values;
    }

    int *getColIdx() const {
        return col_idx;
    }

    void setColIdx(int *colIdx) {
        col_idx = colIdx;
    }

    int *getRowOffsets() const {
        return row_offsets;
    }

    void setRowOffsets(int *rowOffsets) {
        row_offsets = rowOffsets;
    }

    S_Matrix(int row, int col, T *mat) {
        this->row = row;
        this->col = col;
        int cnt = 0;
        for (int i = 0; i < col * row; ++i) {
            if (mat[i] != 0) cnt++;
        }

        values = new T[cnt];
        col_idx = new int[cnt];
        row_offsets = new int[row + 1];

        int j = 0;
        int row_num = 0;
        bool find = false;
        row_offsets[0] = -1;
        for (int i = 0; i < col * row; ++i) {
            if (i / col == row_num + 1) {
                if (!find) row_offsets[row_num] = row_offsets[row_num - 1];
                row_num++;
                find = false;
            }
            if (mat[i] != 0) {
                values[j] = mat[i];
                col_idx[j] = i % col;
                if (!find) {
                    find = true;
                    row_offsets[row_num] = j;
                }
                j++;
            }
        }
        if (!find) row_offsets[row - 1] = cnt;
        row_offsets[row] = cnt;
    }

    ~S_Matrix(){
        delete [] values;
        delete [] col_idx;
        delete [] row_offsets;
    }

    template<typename T1>
    S_Matrix operator*(S_Matrix<T1> &other) {
        T *result = new T[row * other.getCol()];
        for (int i = 0; i < row * other.getCol(); ++i) {
            result[i] = 0;
        }
        for (int i = 0; i < row; ++i) {
            while (row_offsets[i] == row_offsets[i+1]){
                i++;
            }
            if (i>row) break;
            //第i行，第col_idx[j]列，要找other矩阵的第col_idx[j]行，other.col_idx[k]列
            for (int j = row_offsets[i]; j < row_offsets[i+1]; ++j) {
                int begin = other.getRowOffsets()[col_idx[j]];
                int end = other.getRowOffsets()[col_idx[j]+1];
                for (int k = begin; k < end; ++k) {
                    T1 to_mul = other.findElement(col_idx[j],other.getColIdx()[k]);
                    if (to_mul!=0) result[i*other.getCol()+other.getColIdx()[k]] += (values[j]*to_mul);
                }
            }
        }
        return S_Matrix{row,other.getCol(),result};
    }

    T findElement(int row_num, int col_num){
        T res = 0;
        for (int i = row_offsets[row_num]; i < row_offsets[row_num+1]; ++i) {
            if (col_idx[i] == col_num) {
                res = values[i];
                break;
            }
        }
        return res;
    }

    void showInfo() {
        cout << "values: ";
        for (int i = 0; i < row_offsets[row]; ++i) {
            cout << values[i] << " ";
        }
        cout << "\ncol_idx: ";
        for (int i = 0; i < row_offsets[row]; ++i) {
            cout << col_idx[i] << " ";
        }
        cout << "\nrow_offset: ";
        for (int i = 0; i < row + 1; ++i) {
            cout << row_offsets[i]<<" ";
        }
    }
};


#endif //PROJECT_S_MATRIX_HPP
