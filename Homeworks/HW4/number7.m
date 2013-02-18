#!/usr/bin/octave
# Created by Hershal Bhave on 2/16/13
# For M368K HW3, § 8.5 Number 7
# Written in GNU Octave

# Original Function
#f=@(x) x.^2 .* cos(x);

# Domain Extrema
a=-pi;
b=pi;

# Plot Domain
x=a:.01:b;

# Degree of the polynomial
n=2;
m=3;

[A,B] = dtriglsq(f,n,m,a,b);

s=A(1)./2 + A(n+1).*cos(n.*x);
for i=1:n-1
  s+= A(i+1).*cos(i.*x) + B(i+1).*sin(i.*x);
endfor

figure;
hold on;
plot(x,f(x), 'b');
plot(x,s, 'k');
hold off;

legend("Original Function f(x)","Discrete Trig Poly Approx","location","northwest");
title("§ 8.5 Number 7: y vs x of f(x) and its Discrete Least-Squares Trig Poly Approximation")
xlabel("x");
ylabel("y");
