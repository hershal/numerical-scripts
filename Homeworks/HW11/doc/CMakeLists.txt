cmake_minimum_required(VERSION 3.0)
project (hw11-doc)
include (../../UseLATEX.cmake)

set (doc-staging-area "./doc-staging")
set (dependent-files ../scripts/poissonFinDiff.m ../src/program11.cpp ../../libnumerics/linearcd2D.cpp)
file (MAKE_DIRECTORY ${doc-staging-area})

add_custom_target (copy)
add_custom_command (
  TARGET copy
  COMMAND ${CMAKE_COMMAND} -E copy ../../scripts/poissonFinDiff.m ${doc-staging-area}/poissonFinDiff.m
  COMMAND ${CMAKE_COMMAND} -E copy ../../src/program11.cpp ${doc-staging-area}/program11.cpp
  COMMAND ${CMAKE_COMMAND} -E copy ../../../libnumerics/linearcd2D.cpp ${doc-staging-area}/linearcd2D.cpp
)

add_latex_document (
  writeup.tex
  DEPENDS copy
  CONFIGURE writeup.tex
  DEFAULT_PDF
)
