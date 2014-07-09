cmake_minimum_required(VERSION 3.0)
project (hw13-doc)
include (../../UseLATEX.cmake)

set (doc-staging-area "./doc-staging")

add_custom_target (copy)
add_custom_command (
  TARGET copy
   COMMAND ${CMAKE_COMMAND} -E copy ../../include/gnuplot-lua-tikz.sty ./gnuplot-lua-tikz.sty
   COMMAND ${CMAKE_COMMAND} -E copy ../../include/gnuplot-lua-tikz-common.tex ./gnuplot-lua-tikz-common.tex
   COMMAND ${CMAKE_COMMAND} -E copy ../../scripts/contour.pl ${doc-staging-area}/contour.pl
   COMMAND ${CMAKE_COMMAND} -E copy ../../scripts/plot.gp ${doc-staging-area}/plot.gp
   COMMAND ${CMAKE_COMMAND} -E copy ../../scripts/wavediff.m ${doc-staging-area}/wavediff.m
   COMMAND ${CMAKE_COMMAND} -E copy ../../src/program13.cpp ${doc-staging-area}/program13.cpp
   COMMAND ${CMAKE_COMMAND} -E copy ../../../libnumerics/ctrdiff2D.cpp ${doc-staging-area}/ctrdiff2D.cpp
 )

 add_latex_document (
   writeup.tex
   DEPENDS copy
   CONFIGURE writeup.tex
   INPUTS gnuplot-lua-tikz.sty gnuplot-lua-tikz-common.tex
   DEFAULT_PDF
)