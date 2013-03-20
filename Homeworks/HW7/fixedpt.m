#!/usr/bin/octave
# Created by Hershal Bhave on 3/19/13
# For M368K HW3, § 10.1 Number 6
# Written in GNU Octave
#
# Description: Finds the fixed point given an array of x equations
#
# **** WARNING ****
# This is not to be used in a production environment as it does not
# hold up to quality standards. This implies the TODO below.
# **** WARNING ****
#
# TODO: Make this method generalized for functions of n variables
#

function x = fixedpt(x, x1, x2, x3, n)

  i = 1;
  do
      
    # Normal Style
    x(i+1,1) = x1(x(i,1), x(i,2), x(i,3));
    x(i+1,2) = x2(x(i,1), x(i,2), x(i,3));
    x(i+1,3) = x3(x(i,1), x(i,2), x(i,3));

    ## x(i+1,1) = x1(x(i,1), x(i,2));
    ## x(i+1,2) = x2(x(i,1), x(i,2));

    # G-S style
    ## x(i+1,1) = x1(x(i,1), x(i,2), x(i,3));
    ## x(i+1,2) = x2(x(i+1,1), x(i,2), x(i,3));
    ## x(i+1,3) = x3(x(i+1,1), x(i+1,2), x(i,3));

    ## x(i+1,1) = x1(x(i,2));
    ## x(i+1,2) = x2(x(i+1,1));

    i++;
    until(norm(x(i,:)-x(i-1,:), inf) < 10^-5 || i>n)

    if (i>n)
       printf("Maximum Iterations reached. Breaking...\n");
    endif

endfunction
