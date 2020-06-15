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
void testRGBTrans();
void testRGBConvolution();

int main(){
    cout << "OpenCV Version: " << CV_VERSION << endl;
    if(img.empty()){
        cout<<"can not load image \n";
        return -1;
    }

//    testTrans();
//    testConvolution();

    testRGBTrans();
    testRGBConvolution();
}

void testTrans(){
    cvtColor(img, img, CV_BGR2GRAY);
    imshow("Pokemon",img);

    N_Matrix<int> a{img,0};
    N_Matrix<int> b = a.transposition();

    imshow("PokemonT", b.toOpenCVMat(0));

    waitKey(0);
}

void testConvolution(){
    cvtColor(img, img, CV_BGR2GRAY);
    N_Matrix<double> a{img,0};

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
    imshow("MyConvol", b.toOpenCVMat(0));
    imshow("OpecCVConvol", dstImage);

    waitKey(0);
}

void testRGBTrans(){
//    imshow("Pokemon",img);
    N_Matrix<int> a{img,1};
    N_Matrix<int> r,g,b;

    fromRGB<int>(img,r,g,b);

    b = b.transposition();
    g = g.transposition();
    r = r.transposition();

    imshow("Pokemon1",a.toOpenCVMat(1));
    imshow("Pokemon2",toRGB<int>(r,g,b));

    waitKey(0);
}

void testRGBConvolution(){
    N_Matrix<int> a{img,1};
    N_Matrix<int> r,g,b,r1,g1,b1,r2,g2,b2,
                        r3,g3,b3,r4,g4,b4;

    fromRGB<int>(img,r,g,b);

    double kernel1[] = {
            1.0/9,1.0/9,1.0/9,
            1.0/9,1.0/9,1.0/9,
            1.0/9,1.0/9,1.0/9
    };
    double kernel2[] = {
            -1,2,-1,
            0,0,0,
            -1,2,-1
    };    double kernel3[] = {
            -1,-1,-1,
            -1,9,-1,
            -1,-1,-1
    };    double kernel4[] = {
            0,-1,0,
            -1,4,-1,
            0,-1,0
    };
    N_Matrix<double> k1{3,3,kernel1};
    N_Matrix<double> k2{3,3,kernel2};
    N_Matrix<double> k3{3,3,kernel3};
    N_Matrix<double> k4{3,3,kernel4};

    b1 = b.convolution(k1);
    g1 = g.convolution(k1);
    r1 = r.convolution(k1);
    b2 = b.convolution(k2);
    g2 = g.convolution(k2);
    r2 = r.convolution(k2);
    b3 = b.convolution(k3);
    g3 = g.convolution(k3);
    r3 = r.convolution(k3);
    b4 = b.convolution(k4);
    g4 = g.convolution(k4);
    r4 = r.convolution(k4);

    imshow("Pokemon0",toRGB<int>(r,g,b));
    imshow("Pokemon1",toRGB<int>(r1,g1,b1));
    imshow("Pokemon2",toRGB<int>(r2,g2,b2));
    imshow("Pokemon3",toRGB<int>(r3,g3,b3));
    imshow("Pokemon4",toRGB<int>(r4,g4,b4));

    waitKey(0);
}