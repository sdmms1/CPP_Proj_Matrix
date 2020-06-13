#include <iostream>
#include "N_Matrix.hpp"
using namespace std;

int A[] = {1,2,3,
           4,5,6,
           7,8,9};
double B[] = {9,8,7,
              6,5,4,
              3,2,1};
complex<double> C[] = {complex<double>{1,2},complex<double>{7,-5.1},
                       complex<double>{-1,-1},complex<double>{-2.3,1}};

int main() {
    N_Matrix<int> a{3,3,A};
    N_Matrix<double> b{3,3,B};
    N_Matrix<complex<double>> c{2,2,C};

    (a + b).printMatrix();
    (a - b).printMatrix();
    (a.scalar_product(3)).printMatrix();
    (a / 3).printMatrix();
    a.transposition().printMatrix();
    c.conjugation().printMatrix(); //todo
    (a * b).printMatrix();
    (a.element_wise_mult(b)).printMatrix();
    (a.dot_product<double>(b)).printMatrix();   //todo
    (a.crossProduct()).printMatrix();   //todo

    cout << a.getMax() << endl;
    cout << a.getMaxByCol(0) << endl;
    cout << a.getMaxByRow(1) << endl;
    cout << a.getMin() << endl;
    cout << a.getMinByCol(0) << endl;
    cout << a.getMinByRow(1) << endl;
    cout << a.getSum() << endl;
    cout << a.getSumByCol(0) << endl;
    cout << a.getSumByRow(1) << endl;
    cout << a.getAvg() << endl;
    cout << a.getAvgByCol(0) << endl;
    cout << a.getAvgByRow(1) << endl;
}
