#!/usr/bin/octave
# Created by Hershal Bhave on 2/16/13
# For M368K HW3, § 8.6 Number 2b
# Written in GNU Octave

# Original Function
f=@(x) x.*(pi-x);

# Domain Extrema
a=-pi;
b=pi;

# Plot Domain
xx=a:.01:b;
x = (a:(b-a)/(2*n):b)';
x(length(x)) = [];

# Degree of the polynomial
n=4;

[A,B] = fftlsq(f,n,a,b);

s=(A(1) + A(n+1).*cos(n.*x))/2;
for i=1:n-1
  s+= A(i+1).*cos(i.*x) + B(i+1).*sin(i.*x);
endfor

ss=(A(1) + A(n+1).*cos(n.*xx))/2;
for i=1:n-1
  ss+= A(i+1).*cos(i.*xx) + B(i+1).*sin(i.*xx);
endfor


figure;
hold on;
plot(x,f(x), 'b');
plot(x,s, 'o');
plot(xx,ss,'c');
hold off;

legend("Original Function f(x)","FFT Trig Poly Approx","FFT Trig Continuous Poly", "location","northwest");
title("§ 8.6 Number 1c: y vs x of f(x) and its  FFT Least-Squares Trig Poly Approximation (n=4)")
xlabel("x");
ylabel("y");
