#include <iostream>
#include "N_Matrix.hpp"
using namespace std;

int A[] = {1,2,3,
           4,5,6,
           7,8,9};
double B[] = {9,8,7,
              6,5,4,
              3,2,1};
//complex<double> C[] = {complex<double>{1,2},complex<double>{7,-5.1},
//                       complex<double>{-1,-1},complex<double>{-2.3,1}};

double D[] = {5, 6, 1, 2, 3,
           7, 8, 9, 9, 1,
           5, 5, 3, 6, 2,
           1, 5, 3, 6, 8,
           7, 7, 5, 2, 0};

double E[] = {1.2,0.8,0.8,1.2};

double F[] = {1,3.2,
              2.4,-1,
              5.5,-3.1};

double X[] = {1,3,5};
double Y[] = {4,3,2};

int main() {
    N_Matrix<int> a{3,3,A};
    N_Matrix<double> b{3,3,B};
//    N_Matrix<complex<double>> c{2,2,C};
    N_Matrix<double> d{5,5,D};
    N_Matrix<double> e{2,2,E};
    N_Matrix<double> f{3,2,F};
    N_Matrix<double> x{3,1,X};
    N_Matrix<double> y{3,1,Y};

//    a.printMatrix();
//    b.printMatrix();
//    (a + b).printMatrix();
//    (a - b).printMatrix();
//    (a.scalar_product(3)).printMatrix();
//    (a / 3).printMatrix();
//    (b / 3).printMatrix();

//    a.transposition().printMatrix();
//    c.conjugation().printMatrix();
//    (a * b).printMatrix();
//    (a.element_wise_mult(b)).printMatrix();
//    x.dot_product(y).printMatrix();
//    (f.crossProduct()).printMatrix();

//    cout << a.getMax() << endl;
//    cout << a.getMaxByCol(0) << endl;
//    cout << a.getMaxByRow(1) << endl;
//    cout << a.getMin() << endl;
//    cout << a.getMinByCol(0) << endl;
//    cout << a.getMinByRow(1) << endl;
//    cout << a.getSum() << endl;
//    cout << a.getSumByCol(0) << endl;
//    cout << a.getSumByRow(1) << endl;
//    cout << a.getAvg() << endl;
//    cout << a.getAvgByCol(0) << endl;
//    cout << a.getAvgByRow(1) << endl;

//    cout << e.getDet() << endl << endl;
//    d.getInverse().printMatrix();   //todo 1.非法异常   2.自动转换double
//    cout << endl;
//    cout << e.getTrace() << endl << endl;
    a.getEigenValue().printMatrix();
//    cout << endl;
//    e.getEigenvector().printMatrix();
//    a.getEigenvector().printMatrix();

//    d.sliceByCol(1,3).printMatrix();
//    d.sliceByRow(2,3).printMatrix();
//    d.slice(1,2,3,3).printMatrix();
//    f.reShape(2,3);
//    f.printMatrix();

}
