#!/usr/bin/octave
# Created by Hershal Bhave on 3/01/13
# For M368K HW3, § 9.4 Number 2
# Written in GNU Octave
#
# Description: Finds the QR factorization of a matrix
# NOT FINISHED!!
#

function A = qr(A)

B=A;

n=length(A);
G=eye(n,n);
for k=1:n-1

  P=eye(n,n);
  denom=sqrt(A(k,k)^2+A(k+1,k)^2);

  c=A(k,k)/denom;
  s=A(k+1,k)/denom;

  P(k,k)=c;
  P(k+1,k+1)=c;
  P(k,k+1)=s;
  P(k+1,k)=-s;

  A=P*A;
  G*=P';
  P
  A

endfor	 
G
A=A*G


endfunction
