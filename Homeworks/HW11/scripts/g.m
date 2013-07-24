#!/usr/bin/octave
# Created by Hershal Bhave on 04/18/13
# For M368K HW11, § 12.1 Number 2
# Written in GNU Octave
#
# Description: g
#

function z = g(x,y)

  ## For Example 2
  ## if(x==0) 
  ##   z=0;
  ## elseif(x==2)
  ##   z=2.*e.^y;
  ## elseif(y==0)
  ##   z=x;
  ## elseif(y==1)
  ##   z=e.*x;
  ## else
  ##   printf("ERROR IN G\n");
  ## endif

  ## For number 2
  ## if(x==1)
  ##   z=log(y.^2+1);
  ## elseif(x==2)
  ##   z=log(y.^2+4);
  ## elseif(y==0)
  ##   z=2.*log(x);
  ## elseif(y==1)
  ##   z=log(x.^2+1);
  ## else
  ##   printf("ERROR IN G\n");
  ## endif

  ## For number 8
  if(x==0)
    z=y.*(5-y);
  elseif(x==6)
    z=0;
  elseif(y==0)
    z=x.*(6-x);
  elseif(y==5)
    z=0;
  else
    printf("ERROR IN G\n");
  endif

endfunction
