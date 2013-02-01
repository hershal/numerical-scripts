#!/usr/bin/octave

n=5;
A=zeros(n,n);
b=zeros(n,1);
x=zeros(n,1);

w=2.3;

A = diag(ones(n-1,1)*-1,1)+diag(ones(n-1,1)*-1,-1)+diag(2+(1:n)/10,0);
b = [1+(1:n)/10]';

x_sor = sor(A, w, b, x);
x_congrad = congrad(A, eye(n,n), b, x);

fprintf("w: %g\nSOR: %g\nCONGRAD: %g\n",w, size(x_sor,2), size(x_congrad,2));
