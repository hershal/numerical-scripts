#!/usr/bin/octave
# Created by Hershal Bhave on 2/15/13
# For M368K HW3, § 8.5 Number 6
# Written in GNU Octave

f=@(x) number6func(x);

# Plot Domain
x=-pi:.01:pi;

# Degree of the polynomial
n=100;

[a,b] = ctriglsq(f,n,-pi,pi);

s=a(1)./2;
for i=1:n
  s+= a(i+1).*cos(i.*x) + b(i+1).*sin(i.*x);
endfor

figure;
hold on;
plot(x,f(x), 'b');
plot(x,s, 'k');
hold off;

legend("Original Function","Trig Poly Approx","location","northwest");
title("§ 8.5 Number 6: y vs x of f(x) and its Least-Squares Trig Approximation")
xlabel("x");
ylabel("y");
