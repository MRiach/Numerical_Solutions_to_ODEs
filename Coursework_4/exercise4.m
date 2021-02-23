%This produces all plots in exercise 4

%Explicit solution
figure(1)

X = ERK(0,100,[1,0,0],1e-4);
t = 0:1e-4:100;
plot(t,X(:,1))
hold on 
plot(t,X(:,2))
plot(t,X(:,3))
hold off 
title('Numerical Solution to IVP up to t = 100 with 4 stage ERK')
xlabel('$t$','interpreter','latex')
ylabel('$x(t),y(t),z(t)$','interpreter','latex')
legend('$x(t)$','$y(t)$','$z(t)$','Location','southeast','interpreter','latex')

%Implicit solution
figure(2)

X = IRKq4(0,100,[1,0,0],1e-4,5);
t = 0:1e-4:100;
plot(t,X(:,1))
hold on 
plot(t,X(:,2))
plot(t,X(:,3))
hold off 
title('Numerical Solution to IVP up to t = 100 with 4 stage IRK')
xlabel('$t$','interpreter','latex')
ylabel('$x(t),y(t),z(t)$','interpreter','latex')
legend('$x(t)$','$y(t)$','$z(t)$','Location','southeast','interpreter','latex')

