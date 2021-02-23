%Run this to produce plots and global errors for exercise 2
function exercise2()

%This code collects the solutions with different h values using my 
% TaylorSeriesMethod function and then plots them

x1 = TaylorSeriesMethod(0,2,[1,0],0.025);
y1 = 0:0.025:2;
x2 = TaylorSeriesMethod(0,2,[1,0],0.05);
y2 = 0:0.05:2;
x3 = TaylorSeriesMethod(0,2,[1,0],0.1);
y3 = 0:0.1:2;

figure

plot(y1,x1(:,1))

hold on 
plot(y2,x2(:,1))
plot(y3,x3(:,1))
hold off

title('Numerical Solution to IVP 3 using TS(3)')
xlabel('Time, t')
ylabel('Position, x')

legend('h=0.025','h=0.05','h=0.1','Location','northeast')

%This code returns the global error attributed to each h value by making
%use of my TaylorSeriesGlobalError function

h = [0.025,0.05,0.1];
globalerrors = zeros(1,length(h));
for i = 1:length(h)
 globalerrors(i) = TaylorSeriesGlobalError(0,2,[1,0],h(i));
end
X = ['The global error at h = 0.025 is ',num2str(globalerrors(1))];
disp(X)
X = ['The global error at h = 0.05 is ',num2str(globalerrors(2))];
disp(X)
X = ['The global error at h = 0.1 is ',num2str(globalerrors(3))];
disp(X)
end