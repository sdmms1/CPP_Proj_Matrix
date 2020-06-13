//
// Created by sdmms on 2020/6/13.
//
#include <iostream>
#include <opencv2/opencv.hpp>
#include <N_Matrix.hpp>
using namespace std;
using namespace cv;

int main(){
    cout << "OpenCV Version: " << CV_VERSION << endl;
    Mat img = imread("../Pokemon02.png");
    if(img.empty())
    {
        cout<<"can not load image \n";
        return -1;
    }

    cout << img.rows << " " << img.cols << endl;
//    cout << img << endl;
    cvtColor(img, img, CV_BGR2GRAY);
// 求极值 最大、最小值及其位置
    imshow("Pokemon", img);

    int r = img.rows, c = img.cols;
    N_Matrix<int>* m = new N_Matrix<int>{r,c};
    for(int i = 0; i < r; i++){
        for(int j = 0; j < c; j++){
            m->matrix[i*c+j] = img.at<uchar>(i, j);
        }
    }

    cout << img << endl<< endl;
//    m->printMatrix();

    auto n = m->transposition();

    Mat* img2 = new Mat{n.row,n.col,CV_8U, Scalar::all(0)};
    for(int i = 0; i < n.row; i++){
        for(int j = 0; j < n.col; j++){
            img2->at<uchar>(i, j) = n.matrix[i*n.col+j];
        }
    }
    cout << *img2 << endl<< endl;
    imshow("Pokemon2", *img2);
    imshow("Pokemon3",img.t());
    waitKey(0);
    return 0;
};