%This produces all plots in exercise 1


%plot of region of absolute stability for AB(5)
%and return the left hand side of the interval of abs stab

s = linspace(0,2*pi,10000);
r = exp(1i*s);
h = (r.^4-r.^5)./(-1901/720*r.^4+2774/720*r.^3-2616/720*r.^2+1274/720*r-251/720);

figure(1)

plot(h);
title('Locus of AB 5','interpreter','latex')
xlabel('$$Re(\hat{h})$$','interpreter','latex')
ylabel('$$Im(\hat{h})$$','interpreter','latex')
min(real(h));
X = ['The minimum real part of h_hat is ',num2str(min(real(h)))];
disp(X)

%3 plots of numerical solutions vs analytical solutions

X = AB5(0,[1,1],1e-18);
t = 0:1e-18:1e-18*10000;

figure(2)

plot(t,X(:,1))

hold on 
plot(t,ExactSolutionx1(t))
hold off


title('Numerical Solution to IVP, $$u(t)$$, using $$h =\frac{1}{10^{18}} $$','interpreter','latex')
xlabel('Time, t')
ylabel('Position, x')
legend('$x_n$','$x(t_n)$','Location','southwest','interpreter','latex')

X = AB5(0,[1,1],1e-5);
t = 0:1e-5:1e-5*10000;

figure(3)

plot(t,X(:,1))

hold on 
plot(t,ExactSolutionx1(t))
hold off


title('Numerical Solution to IVP, $$u(t)$$, using $$h =\frac{1}{10^{5}} $$','interpreter','latex')
xlabel('Time, t')
ylabel('Position, x')
legend('$x_n$','$x(t_n)$','Location','southwest','interpreter','latex')

X = AB5(0,[1,1],1);
t = 0:1:1*10000;

figure(4)

plot(t,X(:,1))

hold on 
plot(t,ExactSolutionx1(t))
hold off

title('Numerical Solution to IVP, $$u(t)$$, using $$h =1 $$','interpreter','latex')
xlabel('Time, t')
ylabel('Position, x')
legend('$x_n$','$x(t_n)$','Location','southeast','interpreter','latex')