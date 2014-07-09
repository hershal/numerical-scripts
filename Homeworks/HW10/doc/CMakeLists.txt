cmake_minimum_required(VERSION 3.0)
project (hw10-doc)
include (../../UseLATEX.cmake)

set (doc-staging-area "./doc-staging")

add_custom_target (copy)
add_custom_command (
  TARGET copy
  COMMAND ${CMAKE_COMMAND} -E copy ../../scripts/nonlinfindiff.m ${doc-staging-area}/nonlinfindiff.m
  COMMAND ${CMAKE_COMMAND} -E copy ../../scripts/linfindiff.m ${doc-staging-area}/linfindiff.m
  COMMAND ${CMAKE_COMMAND} -E copy ../../src/program10.cpp ${doc-staging-area}/program10.cpp
)

add_latex_document (
  writeup.tex
  DEPENDS copy
  CONFIGURE writeup.tex
  DEFAULT_PDF
)