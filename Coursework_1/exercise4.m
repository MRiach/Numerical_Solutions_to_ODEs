%Run this to produce plots for exercise 4
function exercise4()

%This code collects the solutions for h = 0.01 using my 
% AdamsBashforth3 function and then plots it

X = AdamsBashforth3(0,100,[1,0,0],0.01);

plot3(X(:,1),X(:,2),X(:,3))
title('Numerical Solution to IVP 6 up to t = 100')
xlabel('x(t)')
ylabel('y(t)')
zlabel('z(t)')



end