%This produces all plots in exercise 2


%4 plots of numerical solutions vs analytical solutions
X = LMMq1(0,20,[2,3],1/1498);
t = 0:1/1498:20;

figure(1)

plot(t,X(:,1))

hold on 
plot(t,ExactSolutionx1(t))
hold off
ylim([-2 3])

title('Numerical Solution to IVP, $$u(t)$$, using $$h =\frac{1}{1498} $$','interpreter','latex')
xlabel('Time, t')
ylabel('Position, x')
legend('$x_n$','$x(t_n)$','Location','southwest','interpreter','latex')

figure(2)

plot(t,X(:,2))

hold on 
plot(t,ExactSolutionx2(t))
hold off
ylim([-2 3])

title('Numerical Solution to IVP, $$v(t)$$, using $$h =\frac{1}{1498} $$','interpreter','latex')
xlabel('Time, t')
ylabel('Position, x')
legend('$x_n$','$x(t_n)$','Location','southwest','interpreter','latex')




X = LMMq1(0,20,[2,3],1/1501);
t = 0:1/1501:20;

figure(3)

plot(t,X(:,1))

hold on 
plot(t,ExactSolutionx1(t))
hold off
ylim([-2 3])

title('Numerical Solution to IVP, $$u(t)$$, using $$h =\frac{1}{1501} $$','interpreter','latex')
xlabel('Time, t')
ylabel('Position, x')
legend('$x_n$','$x(t_n)$','Location','southwest','interpreter','latex')

figure(4)

plot(t,X(:,2))

hold on 
plot(t,ExactSolutionx2(t))
hold off
ylim([-2 3])

title('Numerical Solution to IVP, $$v(t)$$, using $$h =\frac{1}{1501} $$','interpreter','latex')
xlabel('Time, t')
ylabel('Position, x')
legend('$x_n$','$x(t_n)$','Location','southwest','interpreter','latex')


%4 plots of global errors for 2 different h values across both solutions,
%u(t) and v(t)

X = LMMGlobalErrorx1(0,20,[2,3],1/1498);
t = 0:1/1498:20;

figure(5)

plot(t,X)

title('Global error of $$u(t)$$ using $$h =\frac{1}{1498} $$','interpreter','latex')
xlabel('Time, t')
ylabel('$|x(t_n)-x_n|$','interpreter','latex')


X = LMMGlobalErrorx2(0,20,[2,3],1/1498);
t = 0:1/1498:20;

figure(8)

plot(t,X)

title('Global error of $$v(t)$$ using $$h =\frac{1}{1498} $$','interpreter','latex')
xlabel('Time, t')
ylabel('$|x(t_n)-x_n|$','interpreter','latex')

X = LMMGlobalErrorx1(0,20,[2,3],1/1501);
t = 0:1/1501:20;

figure(6)

plot(t,X)
ylim([0 0.1])

title('Global error of $$u(t)$$ using $$h =\frac{1}{1501} $$','interpreter','latex')
xlabel('Time, t')
ylabel('$|x(t_n)-x_n|$','interpreter','latex')


X = LMMGlobalErrorx2(0,20,[2,3],1/1501);
t = 0:1/1501:20;

figure(7)

plot(t,X)
ylim([0 0.1])

title('Global error of $$v(t)$$ using $$h =\frac{1}{1501} $$','interpreter','latex')
xlabel('Time, t')
ylabel('$|x(t_n)-x_n|$','interpreter','latex')



