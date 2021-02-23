%This produces all plots in exercise 5

%Numerical solution to BVP (4) as per CW assignment
figure(1)

U = BVPIRKq5(400,100,0.5,0.5,15);
X = 0:0.5:100;
T = 0:0.5:400;
waterfall(T,X,U);
zlim([-0.1 0.1]);
ylim([0 100]);

title('Numerical Solution to BVP with IRK')
xlabel('$t$','interpreter','latex')
ylabel('$x$','interpreter','latex')
zlabel('$u(t,x)$','interpreter','latex')



%Numerical solution to modified BVP 
figure(2)

Z = ModBVPIRKq5(5,100,1,0.1,15);
X = 0:1:100;
T = 0:0.1:5;

waterfall(T,X,Z);
title('Numerical Solution to modified BVP with IRK')
xlabel('$t$','interpreter','latex')
ylabel('$x$','interpreter','latex')
zlabel('$u(t,x)$','interpreter','latex')
zlim([-5 5])

%Analytical solution to modified BVP
figure(3)

Z = ExactSolutionq3(5,100,1,0.1);

waterfall(T,X,Z);
title('Analytical Solution to modified BVP')
xlabel('$t$','interpreter','latex')
ylabel('$x$','interpreter','latex')
zlabel('$\phi (t,x)$','interpreter','latex')
zlim([-5 5])

%Plot of global error vs size of step
figure(4)

deltavalues = [1,0.5,0.25,0.125,0.0625];
globalerrors = zeros(5,1);

for i = 1:length(deltavalues)
    Z = ModBVPIRKq5(5,100,1,deltavalues(i)/10,15);
    P = ExactSolutionq3(5,100,1,deltavalues(i)/10);
    globalerrors(i) = norm(P(:,5)-Z(:,5));
end

plot(deltavalues,globalerrors)
title('Global Error vs Size of Time and Space Step')
xlabel('$\Delta x$','interpreter','latex')
ylabel('$|e_{T}|$','interpreter','latex')
