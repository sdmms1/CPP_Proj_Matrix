//
// Created by Lenovo on 2020/6/8.
//

#ifndef PROJECT_N_MATRIX_HPP
#define PROJECT_N_MATRIX_HPP

#include <iostream>
#include <opencv2/opencv.hpp>
#include "S_Matrix.hpp"
using namespace std;
using namespace cv;

struct SizeException:public exception{
    int type;
    explicit SizeException(int t){type = t;}
    const char * what() const throw() override{
        if (type == 1) return "The size of the two matrices don't support the add operation";
        else if (type == 2) return "The size of the two matrices don't support the subtract operation";
        else if (type == 3) return "The size of the two matrices don't support the element-wise multiply operation";
        else if (type == 5) return "The cross product in n-dimension space should provide n-1 vectors";
        else if(type==4) return "The size of the two matrices don't support the multiply operation";
        else if(type==6) return "The size you provide does not match for reshape operation";
        else if(type==7) return "The matrix size does not support getting determinant operation";
        else if(type==8) return "The matrix size does not support getting trace operation";
        else if(type==9) return "The matrix size does not support cross product operation (need row-column==1";
        else if(type==10) return "Exception! not invertible!";
        else if(type==11) return "Invalid size of matrix!";
        else return "Undefined exception!";
    }
};

struct IndexOutOfRange:public exception{
    int type;
    explicit IndexOutOfRange(int t){type = t;}
    const char* what() const throw(){
        return "Exception! The row (or column) index you provide is out of range!";
    }
};

template <class T>
class N_Matrix {
private:
    int col;
    int row;
    T * matrix;
public:
    int getCol() const {
        return col;
    }

    int getRow() const {
        return row;
    }

    void setSize(int row, int col){
        try{
            if(row <= 0 || col <= 0)
                throw SizeException(11);
        }catch (SizeException& e){
            cerr<<e.what();
            abort();
        }
        this->row = row;
        this->col = col;
        delete[] matrix;
        matrix = new T[row * col];
    }

    T *getMatrix() const {
        return matrix;
    }

    T getValue(int i,int j) const{
        return matrix[i*col+j];
    }

    explicit N_Matrix(int row = 1,int col = 1){
        try{
            if(row <= 0 || col <= 0)
                throw SizeException(11);
        }catch (SizeException& e){
            cerr<<e.what();
            abort();
        }
        this->row = row;
        this->col = col;
        matrix = new T[col*row];
        matrix[0] = 0;
    }
    N_Matrix(int row,int col,const T*Mat):N_Matrix(row,col){
        for (int i = 0; i < row*col; ++i)
            matrix[i] = Mat[i];
    }
    N_Matrix(const N_Matrix& other):N_Matrix(other.getRow(), other.getCol()){
        for(int i = 0; i < row*col; i++){
            matrix[i] = other.matrix[i];
        }
    }
    template <typename T1>
    explicit N_Matrix(const N_Matrix<T1>& other):N_Matrix(other.getRow(), other.getCol()){
        for(int i = 0; i < row*col; i++){
            matrix[i] = other.getMatrix()[i];
        }
    }
    N_Matrix(const S_Matrix<T>& other):N_Matrix(other.getRow(),other.getCol()){
        for (int i = 0; i < row*col; ++i) {
            matrix[i] = 0;
        }

        for (int i = 0; i < row; ++i) {
            while (other.getRowOffsets()[i] == other.getRowOffsets()[i+1]) i++;
            if (i>row) break;
            for (int j = other.getRowOffsets()[i]; j < other.getRowOffsets()[i+1]; ++j) {
                matrix[i*col+other.getColIdx()[j]] = other.getValues()[j];
            }
        }

    }
    explicit N_Matrix(const Mat img, int color):N_Matrix(img.rows, img.cols){
        if(color == 0){
            for(int i = 0; i < row; i++){
                for(int j = 0; j < col; j++){
                    matrix[i*col+j] = img.at<uchar>(i, j);
                }
            }
        }else{
            delete[] matrix;
            col = 3 * col;
            matrix = new T[row * col];
            for(int i = 0; i < row; i++){
                for(int j = 0; j < col / 3; j++){
                    Vec3b pixel = img.at<Vec3b>(i, j);
                    matrix[i*col + 3 * j] = pixel[0];
                    matrix[i*col + 3 * j + 1] = pixel[1];
                    matrix[i*col + 3 * j + 2] = pixel[2];
                }
            }
        }

    }

    ~N_Matrix(){
        delete[] matrix;
    }

    void setMatrix(T*Mat){
        delete[] matrix;
        matrix = new T[row * col];
        for (int i = 0; i < row*col; ++i)
            matrix[i] = Mat[i];
    }
    void modifyVal(int row_num,int col_num,T val){
        matrix[row_num*col+col_num] = val;
    }

    template <typename T1>
    N_Matrix operator + (const N_Matrix<T1>& other) const{
        try{
            if ((this->col != other.getCol()) | (this->row != other.getRow()))
                throw SizeException(1);
        }catch (SizeException& e){
            cerr<<e.what();
            abort();
        }
        N_Matrix result {row,col};
        for (int i = 0; i < col*row; ++i) {
            result.matrix[i] = matrix[i] + other.getMatrix()[i];
        }
        return result;
    }

    template <typename T1>
    N_Matrix operator - (const N_Matrix<T1>& other) const{
        try{
            if ((this->col != other.getCol()) | (this->row != other.getRow()))
                throw SizeException(2);
        }catch (SizeException& e){
            cerr<<e.what();
            abort();
        }
        N_Matrix result {row,col};
        for (int i = 0; i < col*row; ++i) {
            result.matrix[i] = matrix[i] - other.getMatrix()[i];
        }
        return result;
    }

    template <typename T1>
    N_Matrix scalar_product (T1 val) const{
        N_Matrix result {row,col};
        for (int i = 0; i < col*row; ++i) {
            result.matrix[i] = val*matrix[i];
        }
        return result;
    }

    template <typename T1>
    N_Matrix operator / (T1 val){
        N_Matrix result {row,col};
        for (int i = 0; i < col*row; ++i) {
            result.matrix[i] = matrix[i]/val;
        }
        return result;
    }

    N_Matrix& operator = (const N_Matrix& other){
        if(this != &other){
            delete[] matrix;
            row = other.getRow();
            col = other.getCol();
            matrix = new T[row * col];
            for(int i = 0; i < row * col; i++)
                matrix[i] = other.getMatrix()[i];
        }
        return *this;
    }

    bool operator == (const N_Matrix& other){
        if(row != other.getRow() || col != other.getCol())
            return false;
        for(int i = 0; i < row * col; i++){
            if(matrix[i] != other.getMatrix()[i])
                return false;
        }
        return true;
    }

    N_Matrix transposition(){
        N_Matrix result {col,row};
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < col; ++j) {
                result.matrix[j*row+i] = matrix[i*col+j];
            }
        }
        return result;
    }

    N_Matrix<complex<double>> conjugation(){
        N_Matrix<complex<double>> result {row,col};
        for (int i = 0; i < col*row; ++i) {
            result.matrix[i] = conj(matrix[i]);
        }
        return result;
    }

    template <typename T1>
    N_Matrix element_wise_mult(const N_Matrix<T1>& matrix2){
        try{
            if ((row != matrix2.getRow())|(col != matrix2.getCol()))
                throw SizeException(3);
        }catch (SizeException& e){
            cerr<<e.what();
            abort();
        }
        N_Matrix result {row,col};
        for (int i = 0; i < row*col; ++i) {
            result.matrix[i] = matrix[i] * matrix2.getMatrix()[i];
        }
        return result;
    }

    template <typename T1>
    N_Matrix operator * (const N_Matrix<T1>& other) const {
        try{
            if (col != other.getRow()) throw SizeException(4);
        }catch (SizeException& e){
            cerr<<e.what();
            abort();
        }
        N_Matrix result {row,other.getCol()};
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < other.getCol(); ++j) {
                T res = 0;
                for (int k = 0; k < col; ++k) {
                    res = res+ matrix[i*col+k] * other.getMatrix()[j+k*other.getCol()];
                }
                result.matrix[i*other.getCol()+j] = res;
            }
        }
        return result;
    }

    template <typename T1>
    N_Matrix dot_product(const N_Matrix<T1>& other){
        try {
            if ((col != other.getCol()) || (row != other.getRow()))
                throw SizeException(4);
        }catch (SizeException& e){
            cerr<<e.what();
            abort();
        }
        N_Matrix trans = this->transposition();
        return trans * other;
    }

    N_Matrix rotate(){
        N_Matrix result {row,col};
        for (int i = 0; i < row*col; ++i) {
            result.matrix[i] = matrix[col*row-1-i];
        }
        return result;
    }

    template <typename T1>
    N_Matrix convolution(N_Matrix<T1> & core){
        int row_cnt = row+(core.getRow()/2)*2;
        int col_cnt = col+(core.getCol()/2)*2;
        T * augment = new T [row_cnt*col_cnt];

        //make up 0
        for (int i = 0; i < row_cnt; ++i) {
            for (int j = 0; j < col_cnt; ++j) {
                if (i<core.getRow()/2 | i>=core.getRow()/2+row | j<core.getCol()/2 | j>=core.getCol()/2+col) augment[i * col_cnt + j] = 0;
                else augment[i * col_cnt + j] = matrix[(i - core.getRow()/2) * col + (j - core.getCol()/2)];
            }
        }

        N_Matrix result {row,col};
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < col; ++j) {
                T res = 0;
                for (int k = 0; k < core.getRow(); ++k) {
                    for (int l = 0; l < core.getCol(); ++l) {
                        int b_row = i+k;
                        int b_col = j+l;
                        res += core.getMatrix()[k*core.getCol()+l]*augment[(i+k)*col_cnt+j+l];
                    }
                }
                result.matrix[i*col+j] = res;
            }
        }
        return result;
    }

    void printMatrix(){
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < col; ++j) {
                cout<<matrix[i*col+j]<<"\t";
            }
            cout<<endl;
        }
        cout << endl;
    }

    T getMaxByRow(int rowNum){
        try {
            if(rowNum<0||rowNum>row-1) throw IndexOutOfRange(1);
        }catch(IndexOutOfRange& e){
            cerr<<e.what()<<endl;
            abort();
        }
        T max = matrix[col * rowNum];
        for (int i = 0; i < col; i++) {
            if (matrix[col * rowNum + i] > max)
                max = matrix[col * rowNum + i];
        }
        return max;
    }

    T getMaxByCol(int colNum){
        try {
            if(colNum<0||colNum>col-1) throw IndexOutOfRange(1);
        }catch(IndexOutOfRange& e){
            cerr<<e.what()<<endl;
            abort();
        }
        T max = matrix[colNum];
        for (int i = 0; i < row; i++) {
            if (matrix[col * i + colNum] > max)
                max = matrix[col * i + colNum];
        }
        return max;
    }

    T getMax(){
        T max=matrix[0];
        for(int i=0;i<row*col;i++){
            if(matrix[i]>max)
                max=matrix[i];
        }
        return max;
    }

    T getMinByRow(int rowNum){
        try {
            if(rowNum<0||rowNum>row-1) throw IndexOutOfRange(1);
        }catch(IndexOutOfRange& e){
            cerr<<e.what()<<endl;
            abort();
        }
        T max = matrix[col * rowNum];
        for (int i = 0; i < col; i++) {
            if (matrix[col * rowNum + i] > max) {}
            else
                max = matrix[col * rowNum + i];
        }
        return max;
    }

    T getMinByCol(int colNum){
        try {
            if(colNum<0||colNum>col-1) throw IndexOutOfRange(1);
        }catch(IndexOutOfRange& e){
            cerr<<e.what()<<endl;
            abort();
        }
        T max = matrix[colNum];
        for (int i = 0; i < row; i++) {
            if (matrix[col * i + colNum] > max) {}
            else
                max = matrix[col * i + colNum];
        }
        return max;
    }

    T getMin(){
        T max=matrix[0];
        for(int i=0;i<row*col;i++){
            if(matrix[i]>max){}
            else
                max=matrix[i];
        }
        return max;
    }

    T getSumByRow(int rowNum){
        try {
            if(rowNum<0||rowNum>row-1) throw IndexOutOfRange(1);
        }catch(IndexOutOfRange& e){
            cerr<<e.what()<<endl;
            abort();
        }
        T sum = matrix[col * rowNum];
        for (int i = 1; i < col; i++) {
            sum = sum + matrix[col * rowNum + i];
        }
        return sum;
    }

    T getSumByCol(int colNum){
        try {
            if(colNum<0||colNum>col-1) throw IndexOutOfRange(1);
        }
        catch(IndexOutOfRange& e){
            cerr<<e.what()<<endl;
            abort();
        }
        T sum = matrix[colNum];
        for (int i = 1; i < row; i++) {
            sum = sum + matrix[col * i + colNum];
        }
        return sum;
    }

    T getSum(){
        T sum=matrix[0];
        for(int i=1;i<row*col;i++){
            sum=sum+matrix[i];
        }
        return sum;
    }

    T getAvgByRow(int rowNum){
        return getSumByRow(rowNum)/col;
    }

    T getAvgByCol(int colNum){
        return getSumByCol(colNum)/row;
    }

    T getAvg(){
        return getSum()/(row*col);
    }

    void reShape(int newRow,int newCol){
        try {
            if(newRow*newCol!=row*col) throw SizeException(6);
            if(newRow<1||newCol<1) throw IndexOutOfRange(1);
        }catch(SizeException& e){
            cerr<<e.what()<<endl;
            abort();
        }catch(IndexOutOfRange& e){
            cerr<<e.what()<<endl;
            abort();
        }
//        *this = transposition();
        T *newmatrix = new T[newRow * newCol];
        for (int i = 0; i < newRow; i++) {
            for (int j = 0; j < newCol; j++) {
                int index = i * newCol + j + 1;
                int newi = index % col == 0 ? col - 1 : index % col - 1;
                int newj = index % col == 0 ? index / col - 1 : index / col;
                newmatrix[i * newCol + j] = matrix[newj * col + newi];
            }
        }
        row = newRow;
        col = newCol;
        setMatrix(newmatrix);
        delete[] newmatrix;
    }

    N_Matrix sliceByRow(int startRow,int endRow){
        try {
            if(startRow<0||endRow>row-1||endRow-startRow<0)
                throw IndexOutOfRange(1);
        }catch(IndexOutOfRange& e){
            cerr<<e.what()<<endl;
            abort();
        }

        int newRow = endRow - startRow + 1;
        N_Matrix result{newRow, col};
        for (int i = 0; i < newRow; i++) {
            for (int j = 0; j < col; j++) {
                // result.matrix[i * col + j] = matrix[(startRow + i) * col + j];
                result.modifyVal(i,j,matrix[(startRow + i) * col + j]);
            }
        }
        return result;
    }

    N_Matrix sliceByCol(int startCol,int endCol){
        try {
            if(startCol<0||endCol>col-1||endCol-startCol<0)
                throw IndexOutOfRange(1);
        }catch(IndexOutOfRange& e){
            cerr<<e.what()<<endl;
            abort();
        }
        int newCol = endCol - startCol + 1;
        N_Matrix result{row, newCol};
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < newCol; j++) {
                //result.matrix[i * newCol + j] = matrix[i * col + j + startCol];
                result.modifyVal(i,j,matrix[i * col + j + startCol]);
            }
        }
        return result;
    }

    N_Matrix slice(int x1,int y1, int x2 = 0, int y2 = 0){
        int a = x1 < x2? x1 : x2, b = x1 + x2 - a;
        int c = y1 < y2? y1 : y2, d = y1 + y2 - c;
        try {
            if(a < 0 || b > row - 1 || c < 0 || d > col - 1)
                throw IndexOutOfRange(1);
        }catch(IndexOutOfRange& e){
            cerr<<e.what()<<endl;
            abort();
        }

        return sliceByRow(a,b).sliceByCol(c,d);
    }


    T getDet(){
        try {
            if(col!=row) throw SizeException(7);
        }catch(SizeException& e){
            cerr<<e.what()<<endl;
            abort;
        }

        if (col == 1) {
            return matrix[0];
        } else if (col == 2) {
            return matrix[0] * matrix[3] - matrix[1] * matrix[2];
        }
        T det = 0;
        T *tempArr = new T[(col - 1) * (col - 1)];
        for (int t = 0; t < col; t++) {
            for (int i = 0; i < col - 1; i++) {
                for (int j = 0; j < t; j++) {
                    tempArr[(col - 1) * i + j] = matrix[col * (i + 1) + j];
                }
                for (int j = t; j < col - 1; j++) {
                    tempArr[(col - 1) * i + j] = matrix[col * (i + 1) + j + 1];
                }
            }
            N_Matrix temp(col - 1, col - 1, tempArr);
            det = matrix[t] * temp.getDet() * (t % 2 == 0 ? 1 : -1) + det;

        }
        delete[] tempArr;
        return det;
    }

    bool isInvertible(){
        return !(getDet() == 0);
    }

    N_Matrix getInverse(){
        try {
            if(!isInvertible())
                throw SizeException(10);
        }catch(SizeException& e){
            cerr<<e.what()<<endl;
            abort();
        }
        N_Matrix nmatrix(col, col);
        T *tempArr = new T[(col - 1) * (col - 1)];
        for (int i = 0; i < col; i++) {
            for (int j = 0; j < col; j++) {
                for (int m = 0; m < i; m++) {
                    for (int n = 0; n < j; n++) {
                        tempArr[(col - 1) * m + n] = matrix[col * m + n];
                    }
                    for (int n = j; n < col - 1; n++) {
                        tempArr[(col - 1) * m + n] = matrix[col * m + n + 1];
                    }
                }
                for (int m = i; m < col - 1; m++) {
                    for (int n = 0; n < j; n++) {
                        tempArr[(col - 1) * m + n] = matrix[col * (m + 1) + n];
                    }
                    for (int n = j; n < col - 1; n++) {
                        tempArr[(col - 1) * m + n] = matrix[col * (m + 1) + n + 1];
                    }
                }
                N_Matrix temp(col - 1, col - 1, tempArr);
                nmatrix.modifyVal(j, i, temp.getDet() * ((i + j) % 2 == 0 ? 1 : -1));
            }
        }
        delete[] tempArr;
        return nmatrix / getDet();
    }

    N_Matrix QRdecomposition(){
        N_Matrix nmatrixArr[col];
        nmatrixArr[0]=sliceByCol(0,0);
        for(int i=1;i<col;i++){
            nmatrixArr[i]=sliceByCol(i,i);
            N_Matrix temp=nmatrixArr[0].scalar_product((nmatrixArr[0].transposition()*nmatrixArr[i]/((nmatrixArr[0].transposition()*nmatrixArr[0]).getValue(0,0))).getValue(0,0));
            for(int j=1;j<i;j++){
                temp=temp+nmatrixArr[j].scalar_product((nmatrixArr[j].transposition()*nmatrixArr[i]/((nmatrixArr[j].transposition()*nmatrixArr[j]).getValue(0,0))).getValue(0,0));
            }
            nmatrixArr[i]=nmatrixArr[i]-temp;
        }


        for(int i=0;i<col;i++){
            T sum=nmatrixArr[i].getValue(0,0)*nmatrixArr[i].getValue(0,0);
            for(int j=1;j<col;j++){
                sum=sum+nmatrixArr[i].getValue(j,0)*nmatrixArr[i].getValue(j,0);
            }
            sum=sqrt(sum);
            for(int j=0;j<col;j++){
                nmatrixArr[i].modifyVal(j,0,nmatrixArr[i].getValue(j,0)/sum);
            }
        }
        N_Matrix Q(col,col);
        for(int i=0;i<col;i++){
            for(int j=0;j<col;j++){
                //Q.matrix[i*col+j]=nmatrixArr[j].matrix[i];
                Q.modifyVal(i,j,nmatrixArr[j].getValue(i,0));
            }
        }
        return Q;
    }

    T sqrt(T c){
        double factor=1;
        if(c>1.0e-5){
            factor=1;
        }else{
            while(true){
                c=c*10;
                factor=factor*10;
                if(c>1.0e-5){break;}
            }
        }
        T ans = c;
        double error = 1.0e-10;
        while ((ans * ans - c)/c > error || (c - ans * ans)/c > error) {
            ans = (c / ans + ans) / 2;
        }
        return factor==1?ans:ans/sqrt(factor);

    }

    N_Matrix getEigenValue(){
        N_Matrix copy(row,col);
        copy.setMatrix(matrix);
        T det=copy.getDet();
        int offset=0;

        while(!(det-0>0.01)&&((det-0>-0.01))){
            for(int i=0;i<col;i++){
                //copy.matrix[i*col+i]=copy.matrix[i*col+i]+1;
                copy.modifyVal(i,i,copy.getValue(i,i)+1);
            }
            offset++;
            det=copy.getDet();
        }

        int count=0;
        while(count<200){

            N_Matrix Q=copy.QRdecomposition();
            copy=Q.transposition()*copy*Q;
            count++;
        }
        for(int i=0;i<col;i++){
            // copy.matrix[i*col+i]=copy.matrix[i*col+i]+offset*(-1);
            copy.modifyVal(i,i,copy.getValue(i,i)+offset*(-1));
        }

        return copy;
    }

    N_Matrix getEigenvector(){
        N_Matrix I=getIdentity();
        N_Matrix vec(row,col);
        N_Matrix copy(row,col);

        N_Matrix value=getEigenValue();
        for(int t=0;t<col;t++){
            T offset=value.getValue(t,t)+0.1;
            T det=(*this-I.scalar_product(offset)).getDet();

            while(!(det-0>0.01)&&((det-0>-0.01))){
                offset=offset+0.1;
                det=(*this-I.scalar_product(offset)).getDet();
            }
            copy.setMatrix(matrix);
            copy=(copy-I.scalar_product(offset));
            N_Matrix inverse=copy.getInverse();

            N_Matrix V(row,1);
            for(int i=0;i<row;i++){
                //V.matrix[i]=i + 1;
                V.modifyVal(i,0,i+1);
            }
            int count=0;
            while(count<200){
                V=inverse*V;

                T modulu=V.getModulu();
                for(int i=0;i<row;i++){
                    //V.matrix[i]=V.matrix[i]/modulu;
                    V.modifyVal(i,0,V.getValue(i,0)/modulu);
                }
                count++;
            }
            for(int i=0;i<row;i++){
                //vec.matrix[i*col+t]=V.matrix[i];
                vec.modifyVal(i,t,V.getValue(i,0));
            }
        }
        return vec;
    }

    N_Matrix getIdentity(){
        N_Matrix result(row,col);
        for(int i=0;i<row;i++){
            for(int j=0;j<col;j++){
                if(i==j){
                    result.modifyVal(i,j,1);
                }else{
                    result.modifyVal(i,j,0);
                }
            }
        }
        return result;
    }

    T getTrace(){
        try {
            if(col!=row) throw SizeException(8);
        }catch(SizeException& e){
            cerr<<e.what()<<endl;
            abort();
        }
        T trace = matrix[0];
        for (int i = 1; i < col; i++) {
            trace = trace + matrix[i * col + i];
        }
        return trace;
    }

    N_Matrix crossProduct(){
        try {
            if(row-col!=1) throw  SizeException(9);
        }catch(SizeException& e){
            cerr<<e.what()<<endl;
            abort();
        }

        T *vector = new T[row];
        T *tempArr = new T[col * col];
        for (int t = 0; t < row; t++) {
            for (int j = 0; j < col; j++) {
                for (int i = 0; i < t; i++) {
                    tempArr[(col) * i + j] = matrix[col * i + j];
                }
                for (int i = t; i < col; i++) {
                    tempArr[(col) * i + j] = matrix[col * (i + 1) + j];
                }
            }
            N_Matrix temp(col, col, tempArr);
            vector[t] = temp.getDet() * (t % 2 == 0 ? 1 : -1);

        }
        return N_Matrix(row, 1, vector);
    }

    N_Matrix getHessnberg(){
        N_Matrix result=*this;
        N_Matrix x(row,1);
        N_Matrix w(row,1);
        N_Matrix v(row,1);
        N_Matrix H(row,col);
        for(int t=0;t<col-2;t++){
            for(int i=0;i<=t;i++){
                x.modifyVal(i,0,0);
            }

            for(int i=t+1;i<row;i++){
                x.modifyVal(i,0,result.matrix[i*col+t]);
            }
            T length=x.getValue(t+1,0)*x.getValue(t+1,0);
            for(int i=t+2;i<row;i++){
                length=length+x.getValue(i,0)*x.getValue(i,0);
            }
            length=sqrt(length);
            for(int i=0;i<row;i++){
                if(i!=t+1){
                    w.modifyVal(i,0,0);
                }else{
                    w.modifyVal(i,0,length);
                }
            }
            if(x.getValue(0,0)>0) x=x.scalar_product(-1);
            v=w+x;
            if(v.getLength()-0>0.0001) {
                H = v * v.transposition() / ((v.transposition() * v).getValue(0,0));
                H = H.scalar_product(-2);
                for (int i = 0; i < row; i++) {
                    H.modifyVal(i, i, H.matrix[i * col + i] + 1);
                }
                result = H * result * H;
            }
        }
        return result;
    }

    T getModulu(){
        T result=matrix[0]*matrix[0];
        for(int i=1;i<row;i++){
            result=result+matrix[i]*matrix[i];
        }
        return sqrt(result);
    }

    T getAbsolute(T c){
        return c > 0? c : c*(-1);
    }

    Mat toOpenCVMat(int color){
        if(color == 0){
            Mat img{row,col,CV_8U, Scalar::all(0)};
            for(int i = 0; i < row; i++){
                for(int j = 0; j < col; j++){
                    img.at<uchar>(i, j) = matrix[i*col+j] > 255 ? 255:
                                          (matrix[i*col+j] < 0? 0 : (unsigned char)matrix[i*col+j]);
                }
            }
            return img;
        }else{
            Mat img{row,col / 3,CV_8UC3, Scalar::all(0)};
            for(int i = 0; i < row; i++){
                for(int j = 0; j < col / 3; j++){
                    Vec3b pixel;
                    pixel[0] = matrix[i * col + 3 * j] > 255? 255 :
                               (matrix[i * col + 3 * j] < 0? 0 :
                               matrix[i * col + 3 * j]);
                    pixel[1] = matrix[i * col + 3 * j + 1] > 255? 255 :
                               (matrix[i * col + 3 * j + 1] < 0? 0 :
                               matrix[i * col + 3 * j + 1]);
                    pixel[2] = matrix[i * col + 3 * j + 2] > 255? 255 :
                               (matrix[i * col + 3 * j + 2] < 0? 0 :
                               matrix[i * col + 3 * j + 2]);
                    img.at<Vec3b>(i, j) = pixel;
                }
            }
            return img;
        }

    }
};

template <typename T>
void fromRGB(const Mat& img, N_Matrix<T>& r, N_Matrix<T>& g, N_Matrix<T>& b){
    N_Matrix<int> a{img,1};

    int row = a.getRow(), col = a.getCol() / 3;
    r.setSize(row, col);
    g.setSize(row, col);
    b.setSize(row, col);

    for(int i = 0; i < row; i++){
        for(int j = 0; j < col; j++){
            b.modifyVal(i,j,a.getValue(i, 3 * j));
            g.modifyVal(i,j,a.getValue(i, 3 * j + 1));
            r.modifyVal(i,j,a.getValue(i, 3 * j + 2));
        }
    }
}

template <typename T>
Mat toRGB(const N_Matrix<T>& r, const N_Matrix<T>& g, const N_Matrix<T>& b){
    N_Matrix<int> dst{b.getRow(), 3 * b.getCol()};

    for(int i = 0; i < b.getRow(); i++){
        for(int j = 0; j < b.getCol(); j++){
            dst.modifyVal(i, 3 * j, b.getValue(i,j));
            dst.modifyVal(i, 3 * j + 1, g.getValue(i,j));
            dst.modifyVal(i, 3 * j + 2, r.getValue(i,j));
        }
    }
    return dst.toOpenCVMat(1);
}

#endif //PROJECT_N_MATRIX_HPP