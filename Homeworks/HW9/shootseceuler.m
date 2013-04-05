#!/usr/bin/octave
# Created by Hershal Bhave on 04/04/13 
# For M368K HW8, § 11.1 Number 2a and 4a
# Written in GNU Octave
#
# Description: Uses Secant-Euler to approximate the solution to the
# second-order differential equation in f using n steps between a and
# b with values at a and b aa and bb respectively. Guesses for the
# slope are given in t0 and t1 and the number of iterations to guess
# another t in m.
#

function [x,y] = shootseceuler(f,n,a,b,aa,bb,t0,t1,m)

  [x,y] = eulersivp(f,n,a,b,aa,t0);
  ybo = y(n+1);
  tko = t0;
  tk = t1;
  [tmp,y] = eulersivp(f,n,a,b,aa,t1);
 
  for k=3:m
    ttk = tk;
    tyb = y(n+1);
    tk = tk- ((y(n+1)-bb)*(tk-tko))/(y(n+1) - ybo);
    [tmp,y] = eulersivp(f,n,a,b,aa,tk);
  
    tko = ttk;
    ybo = tyb;
  endfor

endfunction
