#!/usr/bin/octave
# Created by Hershal Bhave on 2/08/13
# For M368K HW3, § 8.1 Number 4
# Written in GNU Octave

x=[1 1.1 1.3 1.5 1.9 2.1]';
y=[1.84 1.96 2.21 2.45 2.94 3.18]';

figure;

subplot(2,1,1);
hold on;
plot(x,y,"o");

a=lsq(x,y,2)';
h=x(1):.01:x(length(x));

plot(h,polyval(flipud(a'),h));
hold off;
title("8.1 Number 4: Data Points and LS Poly");
legend("Data Points","LS Poly","location","northwest");
xlabel("x");
ylabel("y");


subplot(2,1,2);
for i=1:length(x)
  e(i)=sum((y(i)-polyval(flipud(a),x(i)))^2);
endfor
plot(x,e,"o-");
title("8.1 Number 4: Error in LS Poly");
legend("Error E","location","northwest");
xlabel("x");
ylabel("Error");
