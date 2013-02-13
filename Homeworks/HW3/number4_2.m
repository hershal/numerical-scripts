#!/usr/bin/octave
# Created by Hershal Bhave on 2/08/13
# For M368K HW3, § 8.2 Number 4
# Written in GNU Octave
#

f=@(x)log(x+2);

a=contlsq(f,2,-1,1);
polyout(flipud(a),'x')
x=0:.01:1;

figure;
hold on;
plot(x,polyval(flipud(a),x),'.k');
plot(x,f(x),'b');
hold off;

title("8.2 Number 4: Original Function and LS Poly");
legend("Continuous L-S Poly","Original Function", "location", "northwest");
xlabel("x");
ylabel("y");
