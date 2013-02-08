#!/usr/bin/octave
# Created by Hershal Bhave on 2/08/13
# For M368K HW3, § 8.1 Number 6d
# Written in GNU Octave

# Description:
# Constructs the least squares approximation in the form be^(ax)
# using the points given in the vectors x and y


function [a, b] = lsexp(x,y)

  n = length(x);

  a=(n*sum(x.*log(y)) - sum(x)*sum(log(y)))/(n*sum(x.^2) - sum(x)^2);
  b=exp((sum(x.^2)*sum(log(y)) - sum(x.*log(y))*sum(x))/(n*sum(x.^2) - sum(x)^2));

endfunction
