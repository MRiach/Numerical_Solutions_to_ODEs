%This produces all plots in exercise 5

%Implicit solution without projection
figure(1)

X = IRKq5noproj(0,50,[-1,0,0.5],1e-4,5);
plot3(X(:,1),X(:,2),X(:,3))
title('Numerical Solution to IVP up to t = 50 with IRK')
xlabel('$x(t)$','interpreter','latex')
ylabel('$y(t)$','interpreter','latex')
zlabel('$z(t)$','interpreter','latex')


%Implicit with projection
figure(2)

X = IRKq5proj(0,50,[-1,0,0.5],1e-4,5);
plot3(X(:,1),X(:,2),X(:,3),'r')
hold on

[X,Y,Z] = sphere;
r = 3;
X2 = X * r;
Y2 = Y * r;
Z2 = Z * r;
surf(X2-1,Y2+3,Z2+0.5)
title('Numerical Solution to IVP up to t = 50 with IRK projected onto sphere')
xlabel('$x(t)$','interpreter','latex')
ylabel('$y(t)$','interpreter','latex')
zlabel('$z(t)$','interpreter','latex')