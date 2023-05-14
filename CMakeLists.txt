cmake_minimum_required(VERSION 3.0)
project(Test)

set(CMAKE_BUILD_TYPE "Release")
#add_definitions("-DENABLE_SSE")
#set(CMAKE_CXX_FLAGS "-std=c++17 ${SSE_FLAGS} -g -O3 -march=native")
set(CMAKE_CXX_FLAGS "-std=c++17 -O3 -fopenmp")
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules) # Add elements to list
#set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)

find_package(OpenCV 4 REQUIRED)
find_package(Sophus REQUIRED)
#find_package(Pangolin REQUIRED) # used to show pointcloud
find_package(CSparse REQUIRED)
find_package(G2O REQUIRED)
find_package(Eigen3 REQUIRED)

include_directories(
        ${OpenCV_INCLUDE_DIRS}
        ${G2O_INCLUDE_DIRS}
        ${Sophus_INCLUDE_DIRS}
        ${EIGEN3_INCLUDE_DIR}
#        "/usr/include/eigen3/"
        ${Pangolin_INCLUDE_DIRS}
        ${CSPARSE_INCLUDE_DIR}
)

# enable drawing plot
add_subdirectory(third_party/matplotlib-cpp)
include_directories(${PROJECT_SOURCE_DIR}/third_party/matplotlib-cpp)

add_executable(test test.cpp)
#target_include_directories(test PUBLIC ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
target_link_libraries(test matplotlib_cpp
        ${OpenCV_LIBS}
        ${Pangolin_LIBRARIES}
        g2o_core g2o_stuff g2o_types_sba g2o_solver_csparse g2o_csparse_extension
        Sophus::Sophus)
