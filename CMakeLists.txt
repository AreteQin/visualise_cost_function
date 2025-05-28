cmake_minimum_required(VERSION 3.0)
project(visualize_cost_function)

set(CMAKE_CXX_FLAGS "-std=c++17")
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules) # Add elements to list

find_package(OpenCV 4 REQUIRED)
find_package(Sophus REQUIRED)
find_package(CSparse REQUIRED)
find_package(G2O REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Glog REQUIRED)

include_directories(
        ${OpenCV_INCLUDE_DIRS}
        ${G2O_INCLUDE_DIRS}
        ${Sophus_INCLUDE_DIRS}
        ${EIGEN3_INCLUDE_DIR}
        ${Pangolin_INCLUDE_DIRS}
        ${CSPARSE_INCLUDE_DIR}
)

# enable drawing plot
add_subdirectory(third_party/matplotplusplus)

add_executable(ab_estimate ab_estimate.cpp)
target_link_libraries(ab_estimate
        ${GLOG_LIBRARIES}
        matplot
        ${OpenCV_LIBS}
        ${Pangolin_LIBRARIES}
        g2o_core g2o_stuff g2o_types_sba g2o_solver_csparse g2o_csparse_extension
        Sophus::Sophus
)

add_executable(a_estimate a_estimate.cpp)
target_link_libraries(a_estimate
        ${GLOG_LIBRARIES}
        matplot
        ${OpenCV_LIBS}
        ${Pangolin_LIBRARIES}
        g2o_core g2o_stuff g2o_types_sba g2o_solver_csparse g2o_csparse_extension
        Sophus::Sophus
)

add_executable(visualize_frames visualize_frames.cpp)
target_link_libraries(visualize_frames
        ${GLOG_LIBRARIES}
        matplot
        ${OpenCV_LIBS}
        ${Pangolin_LIBRARIES}
        g2o_core g2o_stuff g2o_types_sba g2o_solver_csparse g2o_csparse_extension
        Sophus::Sophus
)

add_executable(points_generate points_generate.cpp)
target_link_libraries(points_generate
        ${GLOG_LIBRARIES}
        matplot
        ${OpenCV_LIBS}
        ${Pangolin_LIBRARIES}
        g2o_core g2o_stuff g2o_types_sba g2o_solver_csparse g2o_csparse_extension
        Sophus::Sophus
)
