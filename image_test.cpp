//
// Created by sdmms on 2020/6/13.
//
#include <iostream>
#include "N_Matrix.hpp"
using namespace std;
using namespace cv;

Mat img = imread("../Pokemon02.png");

void testTrans();
void testConvolution();

int main(){
    cout << "OpenCV Version: " << CV_VERSION << endl;
    if(img.empty()){
        cout<<"can not load image \n";
        return -1;
    }

    testTrans();
//    testConvolution();

};

void testTrans(){
    imshow("Pokemon",img);
    cvtColor(img, img, CV_BGR2GRAY);

    N_Matrix<int> a{img};
    N_Matrix<int> b = a.transposition();

    imshow("PokemonG",img);
    imshow("PokemonT", b.toOpenCVMat());

    waitKey(0);
}

void testConvolution(){
    cvtColor(img, img, CV_BGR2GRAY);
    N_Matrix<double> a{img};

    double kernel[] = {
//            1.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/9
            -1,2,-1,0,0,0,-1,2,-1
//            0,-1,0,-1,9,-1,0,-1,0
//            0,-1,0,-1,4,-1,0,-1,0
    };
    N_Matrix<double> k1{3,3,kernel};
    Mat k2 = (Mat_<double>(3,3) <<
                                0,-1,0,-1,4,-1,0,-1,0
                );

    Mat dstImage;
    cv::filter2D(img,dstImage,img.depth(),k2);

    N_Matrix<double> b = a.convolution(k1);

    imshow("OriginPic",img);
    imshow("MyConvol", b.toOpenCVMat());
    imshow("OpecCVConvol", dstImage);

    waitKey(0);
}