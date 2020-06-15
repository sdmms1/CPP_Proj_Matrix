//
// Created by sdmms on 2020/6/15.
//

#include <iostream>
#include "N_Matrix.hpp"
using namespace std;

double A[] = {1,2,0,4,0,6};

int main(){
    S_Matrix<double> a{2,3,A};
    a.showInfo();
}