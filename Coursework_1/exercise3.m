%Run this to produce plots for exercise 3
function exercise3()

%This code collects the solutions with different h values using my 
% AdamsBashforth2 function and then plots them

x1 = AdamsBashforth2(0,100,[0,0],0.025);
y1 = 0:0.025:100;
x2 = AdamsBashforth2(0,100,[0,0],0.05);
y2 = 0:0.05:100;
x3 = AdamsBashforth2(0,100,[0,0],0.1);
y3 = 0:0.1:100;

figure

plot(y1,x1(:,1))
ylim([-0.8 0.8])

hold on 
plot(y2,x2(:,1))
plot(y3,x3(:,1))
hold off

title('Numerical Solution to IVP 5 using AB(2)')
xlabel('Time, t')
ylabel('Position, x')

legend('h=0.025','h=0.05','h=0.1','Location','northwest')


end