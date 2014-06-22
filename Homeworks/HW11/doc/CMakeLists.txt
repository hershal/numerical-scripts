cmake_minimum_required(VERSION 2.8)
project (hw11-doc)
include(../../UseLATEX.cmake)
ADD_LATEX_DOCUMENT(writeup.tex
  INPUTS ../scripts/poissonFinDiff.m ../src/program11.cpp ../../libnumerics/linearcd2D.cpp
  DEFAULT_PDF
)