%This produces all plots in exercise 2

%Numerical solution to BVP (2) as per CW assignment
figure(1)

U = BackEulerq1(10,10,0.01,0.01);
X = 0:0.01:10;
T = 0:0.01:10;
waterfall(X,T,U);

title('Numerical Solution to BVP with Backward Euler Method')
xlabel('$t$','interpreter','latex')
ylabel('$x$','interpreter','latex')
zlabel('$u(t,x)$','interpreter','latex')

