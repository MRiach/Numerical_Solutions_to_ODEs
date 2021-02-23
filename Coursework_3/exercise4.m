%This produces all plots in exercise 4

%plots of loci with using h1 and h2 which are written in terms of r
%using the quadratic formula
figure(1) 

s = linspace(0,2*pi,10000);
r = exp(1i*s);
a = 5/24-5/8*r;
b = 1/12-13/12*r;
c = r.^2-r;

h1 = (-b+sqrt(b.^2-4*a.*c))./(2*a);
h2 = (-b-sqrt(b.^2-4*a.*c))./(2*a);

plot(h1,'b')
hold on
plot(h2,'r')
hold off 

title('Plot of locus where $$|r_1|=1$$','interpreter','latex')
xlabel('$$Re(\hat{h})$$','interpreter','latex')
ylabel('$$Im(\hat{h})$$','interpreter','latex')
legend('$\hat{h}_1$','$\hat{h}_2$','Location','northwest','interpreter','latex')


%plot size of biggest root against h_hat

figure(2)

h = -2.5:0.01:0;
maxr = zeros(length(h),1);
for i = 1:length(h)

    p = [1 -(1+13/12*h(i)+5/8*h(i)^2) h(i)/12+5/24*h(i)^2];
    roots1 = roots(p);
    maxr(i) = max(abs(roots1));
    
end

plot(h,maxr)

title('Plot of $$\hat{h}$$ vs $max(|r_1|,|r_2|)$','interpreter','latex')
xlabel('$$\hat{h}$$','interpreter','latex')
ylabel('$$max(|r_1|,|r_2|)$$','interpreter','latex')


