%This produces all plots in exercise 5

%plot LMM 3 solution
figure(1)

X = LMMq5a(0,100,[1,1,1],0.001);
plot3(X(:,1),X(:,2),X(:,3))
title('Numerical Solution to IVP up to t = 100 with LMM 3-step method')
xlabel('$x(t)$','interpreter','latex')
ylabel('$y(t)$','interpreter','latex')
zlabel('$z(t)$','interpreter','latex')


%plot predictor-corrector method 
figure(2)

X = LMMq5b(0,100,[1,1,1],0.001);
plot3(X(:,1),X(:,2),X(:,3));
title('Numerical Solution to IVP up to t = 100 with predictor-corrector method')
xlabel('$x(t)$','interpreter','latex')
ylabel('$y(t)$','interpreter','latex')
zlabel('$z(t)$','interpreter','latex')
