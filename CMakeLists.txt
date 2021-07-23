# Test on 3.20 
cmake_minimum_required(VERSION 3.20)

# Project name and a few useful settings 
project(
    RiemannProblem
    VERSION 0.1
    DESCRIPTION "2D Riemann problem"
    LANGUAGES Fortran) 
    
# default real size 
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    set(CMAKE_Fortran_FLAGS   "${CMAKE_Fortran_FLAGS} -real-size 64")
endif()

add_subdirectory(src) 