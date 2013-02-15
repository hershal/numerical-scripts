#!/usr/bin/octave
# Created by Hershal Bhave on 2/14/13
# For M368K HW3, § 8.5 Number 2
# Written in GNU Octave

f=@(x) abs(x);


# Plot Domain
x=-pi:.01:pi;

# Degree of the polynomial
n=3;

[a,b] = triglsq(f,n,-pi,pi)

s=a(1)./2 + a(length(a)).*cos((length(a)-1).*x);

for i=2:length(a)-1
  s+= a(i).*cos((i-1).*x) + b(i).*sin((i-1).*x);
endfor

figure;
hold on;
plot(x,f(x), 'b');
plot(x,s, 'k');
hold off;

legend("Original Function","Trig Poly Approx","location","northwest");
title("y vs x of f(x) and its Trig Poly Approximation");
xlabel("x");
ylabel("y");
