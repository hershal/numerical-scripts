cmake_minimum_required(VERSION 3.0)
set (numerics-src matrix.cpp gauss_elim.cpp linearfem.cpp linearcd2D.cpp ctrdiff2D.cpp)
add_library (numerics SHARED ${numerics-src})