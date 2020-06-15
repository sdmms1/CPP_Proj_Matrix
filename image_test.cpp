//
// Created by sdmms on 2020/6/13.
//
#include <iostream>
#include "N_Matrix.hpp"
#include <ostream>
using namespace std;
using namespace cv;

Mat img = imread("../Pokemon02.png");

void testTrans();
void testConvolution();
void testRGB();

int main(){
    cout << "OpenCV Version: " << CV_VERSION << endl;
    if(img.empty()){
        cout<<"can not load image \n";
        return -1;
    }

//    testTrans();
//    testConvolution();

    testRGB();
};

void testRGB(){
    imshow("Image", img);
    for(int i = 0; i < img.rows; i++){
        for(int j = 0; j < img.cols; j++){
            Vec3b pixel = img.at<Vec3b>(i, j);
            pixel[0] = 0;
//            pixel[1] = 0;
//            pixel[2] = 0;
            img.at<Vec3b>(i, j) = pixel;
        }
    }
    imshow("hahaha", img);

    waitKey(0);
}

void testTrans(){
    imshow("Pokemon",img);

    N_Matrix<int> a{img};
    N_Matrix<int> b = a.transposition();

    imshow("PokemonT", b.toOpenCVMat());

    waitKey(0);
}

void testConvolution(){
    N_Matrix<double> a{img};

    double kernel[] = {
//            1.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/9
//            -1,2,-1,0,0,0,-1,2,-1
            -1,-1,-1,
            -1,9,-1,
            -1,-1,-1
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