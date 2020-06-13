#include <iostream>
#include "N_Matrix.hpp"
using namespace std;

int main() {
    int arr[] = {1,2,3,4,5,6,7,8,9};
    N_Matrix<int> a{3,3,arr};

    Mat b = a.toOpenCVMat();
    cout << b << endl;

    N_Matrix<int> c{b};
    c.printMatrix();
}
