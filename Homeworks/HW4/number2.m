#!/usr/bin/octave
# Created by Hershal Bhave on 2/14/13
# For M368K HW3, § 8.5 Number 2
# Written in GNU Octave

# Original Function
f=@(x) x;

# Plot Domain
x=-pi:.01:pi;

# Degree of the polynomial
n=2;

[a,b] = ctriglsq(f,n,-pi,pi);

s=a(1)./2 + a(n+1).*cos(n.*x);
for i=1:n-1
  s+= a(i+1).*cos(i.*x) + b(i+1).*sin(i.*x);
endfor

figure;
hold on;
plot(x,f(x), 'b');
plot(x,s, 'k');
hold off;

legend("Original Function","Trig Poly Approx","location","northwest");
title("§ 8.5 Number 2: y vs x of f(x) and its Least-Squares Trig Poly Approximation");
xlabel("x");
ylabel("y");
