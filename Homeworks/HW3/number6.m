
#x=[1 1.25 1.5 1.75 2]';
#y=[5.1 5.79 6.53 7.45 8.46]';
x=[.2 .3 .6 .9 1.1 1.3 1.4 1.6]
y=[0.050446 0.098426 0.33277 0.72660 1.0972 1.5697 1.8487 2.5015]
h=x(1):.01:x(length(x));

[a1,b1] = lsexp(x,y)
[a2,b2] = lspwr(x,y)
e1=sum((y-b1.*exp(a1.*x)).^2)
e2=sum((y-b2.*x.^a2).^2)

figure;

subplot(2,1,1);
hold on;
plot(x,y,'o');
plot(h,b1*exp(a1*h));
title("8.1 Number 6: Data Points and LS Approximation (b*exp(a*x))");
legend("Data Points","LS Approx","location","northwest");
xlabel("x");
ylabel("y");
hold off;

subplot(2,1,2);
hold on;
plot(x,y,'o');
plot(h,b2.*h.^a2);
title("8.1 Number 6: Data Points and LS Approximation (b*x^a)");
legend("Data Points","LS Approx","location","northwest");
xlabel("x");
ylabel("y");

hold off;

