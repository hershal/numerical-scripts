#!/usr/bin/octave
# Created by Hershal Bhave on 3/25/13
# For M368K HW8, § 10.4 Number 1b
# Written in GNU Octave
#
# Description: Uses the Steepest Descent Method to minimize the
# function f given an initial approximation x and the maximum number
# of steps n
#

function x = steepest(x, f, n)

  tol = 10^-6; k=1;
  gg=@(x)(sum(f(x).^2));

  while k<n
    g(2)=gg(x);
    z=2*jacob(x,f)'*f(x);
    z0=norm(z);
    z=z/z0;

    alpha(2)=0;
    alpha(4)=1;

    g(4)=gg(x-alpha(4).*z);

    while g(4)>g(2)
      alpha(4)=alpha(4)/2;
      g(4)=gg(x-alpha(4).*z);
      if alpha(4)<tol/2
	fprintf("No likely improvement\n");
	k=n;
      endif
    endwhile

    alpha(3)=alpha(4)/2
    g(3)=gg(x-alpha(3).*z)
    
    h1=(g(3)-g(2))/alpha(3)
    h2=(g(4)-g(3))/(alpha(4)-alpha(3))
    h3=(h2-h1)/(alpha(4))

    alpha(1)=0.5*(alpha(3)-h1/h3);

    g(1)=gg(x-alpha(1).*z);
    [gmin,minidx] = min(g);

    x=x-alpha(minidx).*z;

    if abs(gmin-g(2)) < tol
      k=n;
    endif
    k++;
  endwhile

endfunction
