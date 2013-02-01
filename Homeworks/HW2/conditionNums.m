#!/usr/bin/octave
# Created by Hershal Bhave on 1/31/13
# For M368K HW2, § 7.5 Number 8
# Written in GNU Octave

A1=[1 2; 1.00001 2];
b1=[3 3.00001]';

A2=[1 2; 1.000011 2];
b2=[3.00001 3.00003]';

x1 = A1\b1;
x2 = A2\b2;

c1 = norm(A1, inf)*norm(inverse(A1), inf);
c2 = norm(A2, inf)*norm(inverse(A2), inf);

fprintf("A1 solution: %.9g\n", x1);
fprintf("A2 solution: %.9g\n", x2);
fprintf("\n");
fprintf("A1 condition number: %g\nA2 condition number: %g\n", c1, c2);
fprintf("\n");

dA = abs(A1-A2);
db = abs(b1-b2);

lhs = norm(x1-x2)/norm(x1);
rhs = ((condition(A1)*norm(A1))/(norm(A1)-condition(A1)*norm(dA)))*(norm(db)/norm(b1)+norm(dA)/norm(A1));
eq725tf = lhs <= rhs;

if(abs(eq725tf) > 0)
  s = "true";
else
  s = "false";
endif

fprintf("Verifying equation 7.25 in Burden & Faires:\n");
fprintf("%g ?<= %g: %s\n", lhs, rhs, s);
