include(../../UseLATEX.cmake)

ADD_LATEX_DOCUMENT(writeup.tex
  INPUTS ../scripts/linfindiff.m ../scripts/nonlinfindiff.m ../src/program10.cpp
  DEFAULT_PDF
)