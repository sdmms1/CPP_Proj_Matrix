cmake_minimum_required(VERSION 3.16)
project(CPP_Proj_Matrix)

set(CMAKE_CXX_STANDARD 14)

add_executable(CPP_Proj_Matrix main.cpp)

set(OpenCV_INCLUDE_DIRS G:/2020Spring/CS205/LAB/opencv_4_2/mingw_opencv/include/opencv2)
FIND_PACKAGE(OpenCV REQUIRED)
message(STATUS "OpenCV library status:")
message(STATUS "    version: ${OpenCV_VERSION}")
message(STATUS "    libraries: ${OpenCV_LIBS}")
message(STATUS "    include path: ${OpenCV_INCLUDE_DIRS}")
include_directories(
        ${PROJECT_SOURCE_DIR}
        #${OpenCV_INCLUDE_DIRS}
        "G:/2020Spring/CS205/LAB/opencv_4_2/mingw_opencv/include/opencv2"

)
#include_directories(${OpenCV_INCLUDE_DIRS})
add_executable(image_demo image_test.cpp)
target_link_libraries(image_demo ${OpenCV_LIBS})

add_executable(test_demo test1.cpp)
target_link_libraries(test_demo ${OpenCV_LIBS})