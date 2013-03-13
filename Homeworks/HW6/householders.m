#!/usr/bin/octave
# Created by Hershal Bhave on 3/01/13
# For M368K HW3, § 9.4 Number 2
# Written in GNU Octave
#
# Description: Finds the symmetric tridiagonal matrix similar to the
# given symmetric matrix
#

function A = householders(A)

  n = length(A);

  for k=1:n-2

    ## Step 2
    q = 0;
    for j=k+1:n
      q+= (A(j,k,k))^2;
    endfor

    ## Step 3
    if A(k+1,k,k) == 0
      alpha = -q^(1/2);
    else 
      alpha = -((q^(1/2)*A(k+1,k,k))/(abs(A(k+1,k,k))));
    endif

    ## Step 4
    rsq = alpha^2 - alpha*A(k+1,k,k);

    ## Step 5
    v(k) = 0; 
    v(k+1) = A(k+1,k,k) - alpha;

    for j=k+2:n
	v(j) = A(j,k,k);
    endfor

    ## Step 6
    for j=k:n
      ## tmp = 0;
      ## for i=k+1:n
      ## 	tmp += A(j,i,k)*v(i);
      ## endfor
      ## u(j) = (1/rsq)*tmp;
      w(j) = v(j)/sqrt(2*rsq);
    endfor


    P(:,:,k) = eye(size(A,2)) - 2*w'*w;
    A(:,:,k+1) = P(:,:,k)*A(:,:,k)*P(:,:,k);
    
    
    ## ## Step 7
    ## product = 0;
    ## for i=k+1:n
    ##   product += v(i)*u(i);
    ## endfor
    ## 
    ## ## Step 8
    ## for j=k:n
    ##   z(j) = u(j) - (product/(2*rsq))*v(j);
    ## endfor
    ## 
    ## ## Step 9
    ## for l=k+1:n-1
    ##   ## Step 10
    ##   for j=l+1:n
    ## 	A(j,l,k+1) = A(j,l,k) - v(l)*z(j) - v(j)*z(l);
    ## 	A(l,j,k+1) = A(j,l,k+1);
    ##   endfor
    ##   ## Step 11
    ##   A(l,l,k+1) = A(l,l,k) - 2*v(l)*z(l);
    ## endfor
    ## 
    ## ## Step 12
    ## A(n,n,k+1) = A(n,n,k) - 2*v(n)*z(n);
    ## 
    ## ## Step 13
    ## for j=k+2:n
    ##   A(k,j,k+1) = A(j,k,k+1) = 0;
    ## endfor
    ## 
    ## ## Step 14
    ## A(k+1,k,k+1) = A(k+1,k,k) - v(k+1)*z(k);
    ## A(k,k+1,k+1) = A(k+1,k,k);

  endfor

endfunction
