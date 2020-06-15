#include <iostream>
#include <sys\timeb.h>
#include "N_Matrix.hpp"
#define ll long long
using namespace std;

long long getSystemTime(){
    timeb t{};
    ftime(&t);
    return t.time * 1000 + t.millitm;
}

int main() {
    ifstream input;
    input.open("matrix1.txt");
    int mat1 [300][300];
    for (int i = 0; i < 300; ++i) {
        for (int j = 0; j < 300; ++j) {
            input>>mat1[i][j];
        }
    }
    S_Matrix<int> smat1 {300,300,*mat1};
    input.close();
    input.open("matrix2.txt");
    for (int i = 0; i < 300; ++i) {
        for (int j = 0; j < 300; ++j) {
            input>>mat1[i][j];
        }
    }
    S_Matrix<int> smat2 {300,300,*mat1};
    smat1.showInfo();
    cout<<"\n";
    smat2.showInfo();
    cout<<"\n";
    ll start = getSystemTime();
    S_Matrix<int> result = smat1*smat2;
    ll end = getSystemTime();
    result.showInfo();
    cout<<endl;
    cout<<"The sparse matrix multiplication takes "<<end-start<<" ms.";

    N_Matrix<int> nmat1 {smat1};
    N_Matrix<int> nmat2 {smat2};
    start = getSystemTime();
    N_Matrix<int> result2 = nmat1*nmat2;
    end = getSystemTime();
    cout<<endl;
    cout<<"The normal matrix multiplication takes "<<end-start<<" ms.";
    return 0;
}
