cmake_minimum_required(VERSION 2.8)
set (numerics-src matrix.cpp gauss_elim.cpp)
add_library (numerics SHARED ${numerics-src})