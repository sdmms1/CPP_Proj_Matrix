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

double D[] = {1.2,0.8,
              0.8,1.2};

double E[] = {2,5,2,
              3,2,8,
              1,2,3};

double F[] = {8, 6, 1, 5,
              3, 8, 1, 6,
              1, 2, 6, 3,
              2, 5, 6, 4};

double X[] = {1,3,5};
double Y[] = {4,3,2};
double Z[] = {1,4,
              3,3,
              5,2};

template <class T>
void showEigen(N_Matrix<T> a){
    a.printMatrix();
    cout << "Det:" << a.getDet() <<
        ", Trace:" << a.getTrace() << endl << endl;
    cout << "Inverse:" << endl;
    a.getInverse().printMatrix();
    cout << "Eigenvalue:" << endl;
    a.getEigenValue().printMatrix();
    cout << "Eigenvector:" << endl;
    a.getEigenvector().printMatrix();
}

int main() {
    N_Matrix<int> a{3,3,A};
    N_Matrix<double> b{3,3,B};
    N_Matrix<complex<double>> c{2,2,C};
    N_Matrix<double> d{2,2,D};
    N_Matrix<double> e{3,3,E};
    N_Matrix<double> f{4,4,F};
    N_Matrix<double> x{3,1,X};
    N_Matrix<double> y{3,1,Y};
    N_Matrix<double> z{3,2,Y};

//    a.printMatrix();
//    b.printMatrix();
//    (a + b).printMatrix();
//    (a - b).printMatrix();
//    (a.scalar_product(3)).printMatrix();
//    (a / 3).printMatrix();
//    (b / 3).printMatrix();

    a.printMatrix();
    b.printMatrix();
    a.transposition().printMatrix();
    (a * b).printMatrix();
    (a.element_wise_mult(b)).printMatrix();

    c.printMatrix();
    c.conjugation().printMatrix();

    x.printMatrix();
    y.printMatrix();
    x.dot_product(y).printMatrix();

    z.printMatrix();
    (z.crossProduct()).printMatrix();

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


//    showEigen(d);
//    showEigen(e);
//    showEigen(f);

//    f.printMatrix();
//    f.sliceByCol(1,3).printMatrix();
//    f.sliceByRow(2,3).printMatrix();
//    f.slice(1,2,3,3).printMatrix();
//    f.reShape(1,16);
//    f.printMatrix();
//    f.reShape(8,2);
//    f.printMatrix();

}
