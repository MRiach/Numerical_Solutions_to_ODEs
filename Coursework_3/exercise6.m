%This produces all plots in exercise 6

%plot Implicit Method with Newton's method solution
figure(1)

%Take M=5
X = LMMq6(0,100,[1,1,1],0.001,5);
plot3(X(:,1),X(:,2),X(:,3));
title('Numerical Solution to IVP up to t = 100 with Implicit 3-step method')
xlabel('$x(t)$','interpreter','latex')
ylabel('$y(t)$','interpreter','latex')
zlabel('$z(t)$','interpreter','latex')
