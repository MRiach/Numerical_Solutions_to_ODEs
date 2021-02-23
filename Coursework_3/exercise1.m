%This produces all plots in exercise 1


%plot first and second root of characteristic polynomial to see which 
%alpha_2 values satisfy the constraints of convergence
figure(1)

%first root
x = linspace(-5,5,500);
y1 = abs(-1-x+sqrt((x-3).*(x+1)))/2;
y2 = 0*x+1;
plot(x,y1)
hold on
plot(x,y2)
hold off
title('$$|r_1|$$ vs $$\alpha_2$$ with $$|r_1|=1$$','interpreter','latex')
xlabel('$$\alpha_2$$','interpreter','latex')
ylabel('$$|r_1|$$','interpreter','latex')

figure(2) 

%second root
x = linspace(-5,5,500);
y1 = abs(-1-x-sqrt((x-3).*(x+1)))/2;
y2 = 0*x+1;
plot(x,y1)
hold on
plot(x,y2)
hold off
title('$$|r_2|$$ vs $$\alpha_2$$ with $$|r_2|=1$$','interpreter','latex')
xlabel('$$\alpha_2$$','interpreter','latex')
ylabel('$$|r_2|$$','interpreter','latex')


%plots of loci alpha_2 = -1 and -0.5 
figure(3) 

x = linspace(0,2*pi,1000);
h = (exp(3i*x)+(-1)*exp(2i*x)-(1+(-1)))/(3+2*(-1));
plot(h)
title('Locus produced by $$\alpha_2 = -1$$','interpreter','latex')
xlabel('$$Re(\hat{h})$$','interpreter','latex')
ylabel('$$Im(\hat{h})$$','interpreter','latex')

figure(4) 

x = linspace(0,2*pi,1000);
h = (exp(3i*x)+(-0.5)*exp(2i*x)-(1+(-0.5)))/(3+2*(-0.5));
plot(h)
title('Locus produced by $$\alpha_2 = -0.5$$','interpreter','latex')
xlabel('$$Re(\hat{h})$$','interpreter','latex')
ylabel('$$Im(\hat{h})$$','interpreter','latex')


%Code to identify the alpha value which produces the largest interval
Alphas = linspace(-1.5,-0.01,10000);
x = linspace(0,2*pi,2000);
Realmin = zeros(10000);
for i = 1:length(Alphas)
    h = (exp(3i*x)+Alphas(i)*exp(2i*x)-(1+Alphas(i)))/(3+2*Alphas(i));
    for j = 2:length(h)
        %extract first point at which imaginary component becomes 0 again
        %and then take the corresponding real component
        if imag(h(j))<=0
            Realmin(i) = real(h(j-1));
            break
        end
    end
end

%plot of alpha value vs size of interval
figure(5) 

plot(Alphas,abs(Realmin))
title('Length of interval of absolute stability vs $$\alpha_2$$','interpreter','latex')
xlabel('$$\alpha_2$$','interpreter','latex')
ylabel('Length of interval of absolute stability','interpreter','latex')
