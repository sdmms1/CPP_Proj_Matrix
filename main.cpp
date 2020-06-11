#include <iostream>
#include <complex>
#include "N_Matrix.hpp"
#include "S_Matrix.hpp"
using namespace std;

//used to test the two classes
int main() {
    int one[3][2] {1,2,3,4,5,6};
    complex<int> three[3][2] {{{1,1},{2,2}},{{3,3},{4,4}},{{5,5},{6,6}}};
    double two[3][2] {1.5,2.5,3.5,4.5,5.5,6.5};
    int four[2][3] = {1,2,3,4,5,6};
    int five[4][4] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
    int seven[3][3] = {1,2,1,0,0,0,-1,-2,-1};
    N_Matrix<int> a {3,2,*one};
    N_Matrix<double> b {3,2,*two};
    N_Matrix<complex<int>> c {3,2,*three};
    N_Matrix<int> d {2,3,*four};
    N_Matrix<int> h = {4,4,*five};
    N_Matrix<int> e = {3,3,*seven};
    N_Matrix<int> k = h.convolution(e);
    k.printMatrix();

    int s_one[4][4] = {0,0,0,0,0,0,0,0,1,0,7,0,0,0,0,0};
    int s_two[3][3] = {1,7,0,2,0,0,0,0,0};
    int s_three[3][3] = {0,0,3,0,6,0,5,0,0};
    S_Matrix<int> f {3,3,*s_two};
    S_Matrix<int> g {3,3,*s_three};
    S_Matrix<int> i = f*g;
    N_Matrix<int> j {i};
    j.printMatrix();
    return 0;
}
