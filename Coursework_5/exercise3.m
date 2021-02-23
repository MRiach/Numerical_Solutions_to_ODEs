%This produces all plots in exercise 3

%Numerical solution to modified BVP 
figure(1)

Z = BackEulerq3(10,10,0.1,0.1);
X = 0:0.1:10;
T = 0:0.1:10;

waterfall(X,T,Z);
title('Numerical Solution to modified BVP with Backward Euler Method')
xlabel('$t$','interpreter','latex')
ylabel('$x$','interpreter','latex')
zlabel('$u(t,x)$','interpreter','latex')
zlim([-10 10])

%Analytical solution to modified BVP
figure(2)

Z = ExactSolutionq3(10,10,0.1,0.1);

waterfall(X,T,Z);
title('Analytical Solution to modified BVP')
xlabel('$t$','interpreter','latex')
ylabel('$x$','interpreter','latex')
zlabel('$\phi (t,x)$','interpreter','latex')
zlim([-10 10])

%Plot of global error vs size of step
figure(3)

deltavalues = [0.4,0.2,0.1,0.05,0.025];
globalerrors = zeros(5,1);

for i = 1:length(deltavalues)
    Z = BackEulerq3(10,10,deltavalues(i),deltavalues(i));
    P = ExactSolutionq3(10,10,deltavalues(i),deltavalues(i));
    globalerrors(i) = norm(P(:,end)-Z(:,end));
end

plot(deltavalues,globalerrors)
title('Global Error vs Size of Time and Space Step')
xlabel('$\Delta x = \Delta t$','interpreter','latex')
ylabel('$|e_{T}|$','interpreter','latex')

