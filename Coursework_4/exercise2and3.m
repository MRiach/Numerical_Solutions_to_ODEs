%This produces all plots in exercise 2 & 3

%plot to work out interval of abs stability of implicit method
figure(1)

X = linspace(-3,0,10000);
Y = 1+1/8*(7*X+3*X.^2+X.^3+(X-2/3*X.^3-1/3*X.^4)./(1-X));
Y1 = 0*X-1;
plot(X,Y)
title('Plot of stability function for real values of $\hat{h}$','interpreter','latex')
ylabel('$$R(\hat{h})$$','interpreter','latex')
xlabel('$$\hat{h}$$','interpreter','latex')
hold on
plot(X,Y1)
hold off
legend('$R(\hat{h})$','$y=-1$','Location','southeast','interpreter','latex')

%plot to work out region of abs stability of implicit method
%This numerically computes the stability function at a grid of points in
%the complex plane and outputs the values for which the modulus is equal to
%1.
figure(2)

X = linspace(-3,1,1000);
Z=[];
for Y = linspace(-3,3,1000);
   Z = [Z complex(X,Y*ones(1,1000))];
end
Z = Z(1-1e-2<abs(1+1/8*(7*Z+3*Z.^2+Z.^3+(Z-2/3*Z.^3-1/3*Z.^4)./(1-Z))) &...
            abs(1+1/8*(7*Z+3*Z.^2+Z.^3+(Z-2/3*Z.^3-1/3*Z.^4)./(1-Z)))<1+1e-2);
plot(Z,'.');
title('Region of absolute stability','interpreter','latex')
ylabel('$$Im(\hat{h})$$','interpreter','latex')
xlabel('$$Re(\hat{h})$$','interpreter','latex')


%plot to work out interval of abs stability of explicit method
figure(3)

X = linspace(-3,0,10000);
Y = 1+X+1/2*X.^2+1/6*X.^3+1/24*X.^4;
Y1 = 0*X+1;
plot(X,Y)
title('Plot of stability function for real values of $\hat{h}$','interpreter','latex')
ylabel('$$R(\hat{h})$$','interpreter','latex')
xlabel('$$\hat{h}$$','interpreter','latex')
hold on
plot(X,Y1)
hold off
legend('$R(\hat{h})$','$y=1$','Location','southeast','interpreter','latex')

%plot to work out region of abs stability of explicit method
%This numerically computes the stability function at a grid of points in
%the complex plane and outputs the values for which the modulus is equal to
%1.
figure(4)

X = linspace(-3,1,1000);
Z=[];
for Y = linspace(-3,3,1000);
   Z = [Z complex(X,Y*ones(1,1000))];
end
Z = Z(1-1e-2<abs(1+Z+1/2*Z.^2+1/6*Z.^3+1/24*Z.^4) &...
            abs(1+Z+1/2*Z.^2+1/6*Z.^3+1/24*Z.^4)<1+1e-2);
plot(Z,'.');
title('Region of absolute stability','interpreter','latex')
ylabel('$$Im(\hat{h})$$','interpreter','latex')
xlabel('$$Re(\hat{h})$$','interpreter','latex')