%This function plots all figures in exercise 3 
function exercise3()

X = Explicit(0,20,[2,3],0.0001);
Y = Implicit(0,20,[2,3],0.0001);
Z = 0:0.0001:20;

figure

plot(Z,X(:,1))

hold on 
plot(Z,X(:,2))
hold off

title('Numerical Solution to IVP 1 using Explicit LMM')
xlabel('Time, t')
ylabel('Position, x')

legend('u(t)','v(t)','Location','northeast')

figure

plot(Z,Y(:,1))
ylim([-1.5 3])

hold on 
plot(Z,Y(:,2))
hold off

title('Numerical Solution to IVP 1 using Implicit LMM')
xlabel('Time, t')
ylabel('Position, x')

legend('u(t)','v(t)','Location','northeast')