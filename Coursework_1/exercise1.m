%Run this to produce plots and global errors for exercise 1 
function exercise1()

%This code collects the solutions with different h values using my 
% EulerMethod function and then plots them

x1 = EulerMethod(0,3,[0,1],0.025);
y1 = 0:0.025:3;
x2 = EulerMethod(0,3,[0,1],0.05);
y2 = 0:0.05:3;
x3 = EulerMethod(0,3,[0,1],0.1);
y3 = 0:0.1:3;

figure

plot(y1,x1(:,1))

hold on 
plot(y2,x2(:,1))
plot(y3,x3(:,1))
hold off

title('Numerical Solution to IVP 1 using Euler Method')
xlabel('Time, t')
ylabel('Position, x')

legend('h=0.025','h=0.05','h=0.1','Location','northwest')

%This code returns the global error attributed to each h value by making
%use of my EulerGlobalError function
h = [0.025,0.05,0.1];
globalerrors = zeros(1,length(h));
for i = 1:length(h)
 globalerrors(i) = EulerGlobalError(0,3,[0,1],h(i));
end
X = ['The global error at h = 0.025 is ',num2str(globalerrors(1))];
disp(X)
X = ['The global error at h = 0.05 is ',num2str(globalerrors(2))];
disp(X)
X = ['The global error at h = 0.1 is ',num2str(globalerrors(3))];
disp(X)

%This code plots the global error vs the number of time steps which is
%varied from 1 to 5000

timesteps = 1:1:5000;
globalerrors1 = zeros(1,5000);
for i = 1:5000
 globalerrors1(i) = EulerGlobalError(0,3,[0,1],3/i);
end
figure
plot(timesteps,globalerrors1)
ylim([0.0 0.0002])
xlim([0 5000])
hold on 
x=1:5000;
y=0.0001;
plot(x,y*ones(size(x)))
x=3393;
y=0:0.0001:0.0002;
plot(x*ones(size(y)),y)
hold off
title('Global Error vs Number of Time Steps')
xlabel('Number of Time Steps')
ylabel('Global Error')
legend('Global Error','Global Error = 0.0001','Time Steps = 3393','Location','southwest')
end







%