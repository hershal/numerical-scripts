#!/usr/bin/octave
# Created by Hershal Bhave on 1/31/13
# For M368K HW2, § 7.5 Number 8
# Written in GNU Octave

function con = condition(A)

  con = norm(A, inf)*norm(inverse(A), inf);

endfunction

