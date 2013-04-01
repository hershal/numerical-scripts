#!/usr/bin/octave
# Created by Hershal Bhave on 3/25/13
# For M368K HW8, § 10.2 Number 2
# Written in GNU Octave
#
# Description: Derives the function at some tolerance away for
# accurate results
#

function dfdx = derive(f, x)

tol = 10^-6;
## t = [x-tol x x+tol];

for i=1:length(f)
  dfdx(i) = diff(f(x))/diff([x(i)-tol; x(i); x(i)+tol;]);
endfor

endfunction
