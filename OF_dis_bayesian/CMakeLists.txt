cmake_minimum_required (VERSION 2.8) 
project (FastPatchOF) 

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wno-unknown-pragmas -Wall -std=c++11 -msse4")  #-Wall
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wno-unknown-pragmas -Wall -msse4")  #-Wall


FIND_PACKAGE(OpenCV REQUIRED)
FIND_PACKAGE(Eigen3 REQUIRED) 
include_directories(${EIGEN3_INCLUDE_DIR})


set(CODEFILES run_dense.cpp oflow.cpp patch.cpp patchgrid.cpp FDF1.0.1/image.c FDF1.0.1/opticalflow_aux.c FDF1.0.1/solver.c)


add_executable (run_DE_INT ${CODEFILES})
set_target_properties (run_DE_INT PROPERTIES COMPILE_DEFINITIONS "SELECTMODE=2")
set_property(TARGET run_DE_INT APPEND PROPERTY COMPILE_DEFINITIONS "SELECTCHANNEL=1")
TARGET_LINK_LIBRARIES(run_DE_INT ${OpenCV_LIBS})



