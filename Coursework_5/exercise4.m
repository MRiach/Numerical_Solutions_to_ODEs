%This produces all plots in exercise 4

%Numerical solution to BVP 
figure(1)

Z = BVPIRKq4(10,10,0.1,0.01,1,5);
X = 0:0.1:10;
T = 0:0.01:10;

waterfall(T,X,Z);
title('Numerical Solution to BVP with 4 stage IRK and $v=1$','interpreter','latex')
xlabel('$t$','interpreter','latex')
ylabel('$x$','interpreter','latex')
zlabel('$u(t,x)$','interpreter','latex')


figure(2)

Z = BVPIRKq4(10,10,0.1,0.01,0,5);
X = 0:0.1:10;
T = 0:0.01:10;

waterfall(T,X,Z);
title('Numerical Solution to BVP with 4 stage IRK and $v=0$','interpreter','latex')
xlabel('$t$','interpreter','latex')
ylabel('$x$','interpreter','latex')
zlabel('$u(t,x)$','interpreter','latex')


figure(3)


Z = BVPIRKq4(10,10,0.1,0.01,-1,5);
X = 0:0.1:10;
T = 0:0.01:10;

waterfall(T,X,Z);
title('Numerical Solution to BVP with 4 stage IRK and $v=-1$','interpreter','latex')
xlabel('$t$','interpreter','latex')
ylabel('$x$','interpreter','latex')
zlabel('$u(t,x)$','interpreter','latex')


%Global errors

figure(4)

deltavalues = [0.4,0.2,0.1];
globalerrors = zeros(3,1);

for i = 1:length(deltavalues)
    Z = ModBVPIRKq4(10,10,deltavalues(i),deltavalues(i)/200,-1,5);
    P = ExactSolutionq4(10,10,deltavalues(i),deltavalues(i)/200);
    globalerrors(i) = norm(P(:,end)-Z(:,end));
end

plot(deltavalues,globalerrors)
title('Global Error vs Size of Time and Space Step, $v=-1$','interpreter','latex')
xlabel('$\Delta x$','interpreter','latex')
ylabel('$|e_{T}|$','interpreter','latex')

figure(5)

deltavalues = [0.4,0.2,0.1];
globalerrors = zeros(3,1);

for i = 1:length(deltavalues)
    Z = ModBVPIRKq4(10,10,deltavalues(i),deltavalues(i)/200,0,5);
    P = ExactSolutionq4(10,10,deltavalues(i),deltavalues(i)/200);
    globalerrors(i) = norm(P(:,end)-Z(:,end));
end

plot(deltavalues,globalerrors)
title('Global Error vs Size of Time and Space Step, $v=0$','interpreter','latex')
xlabel('$\Delta x$','interpreter','latex')
ylabel('$|e_{T}|$','interpreter','latex')

figure(6)

deltavalues = [0.4,0.2,0.1];
globalerrors = zeros(3,1);

for i = 1:length(deltavalues)
    Z = ModBVPIRKq4(10,10,deltavalues(i),deltavalues(i)/200,1,5);
    P = ExactSolutionq4(10,10,deltavalues(i),deltavalues(i)/200);
    globalerrors(i) = norm(P(:,end)-Z(:,end));
end

plot(deltavalues,globalerrors)
title('Global Error vs Size of Time and Space Step, $v=1$','interpreter','latex')
xlabel('$\Delta x$','interpreter','latex')
ylabel('$|e_{T}|$','interpreter','latex')


%Numerical solution to modified BVP, v = 1
figure(7)

Z = ModBVPIRKq4(10,10,0.1,0.1/200,1,5);
X = 0:0.1:10;
T = 0:0.1/200:10;

waterfall(T,X,Z);
title('Numerical Solution to modified BVP with IRK, $v=1$','interpreter','latex')
xlabel('$t$','interpreter','latex')
ylabel('$x$','interpreter','latex')
zlabel('$u(t,x)$','interpreter','latex')
zlim([-10 10])

%Analytical solution to modified BVP, v = 1
figure(8)

Z = ExactSolutionq4(10,10,0.1,0.1/200);
waterfall(T,X,Z);
title('Analytical Solution to modified BVP, $v=1$','interpreter','latex')
xlabel('$t$','interpreter','latex')
ylabel('$x$','interpreter','latex')
zlabel('$\phi (t,x)$','interpreter','latex')