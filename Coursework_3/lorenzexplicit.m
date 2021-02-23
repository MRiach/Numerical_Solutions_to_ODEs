function w=lorenzexplicit(h)
w=zeros(length(h),3);
t=0:h:100;
sigma=10;
beta=8/3;
rho=28;
w(1,:)=[1,1,1];
w(2,:)=w(1,:)+h*[sigma*(w(1,2)-w(1,1)),w(1,1)*(rho-w(1,3))-w(1,2),w(1,1)*w(1,2)-beta*w(1,3)];
w(3,:)=w(2,:)+h*[sigma*(w(2,2)-w(2,1)),w(2,1)*(rho-w(2,3))-w(2,2),w(2,1)*w(2,2)-beta*w(2,3)];
for n=3:length(t)-1
    w(n+1,1)=0*w(n,1)+1*w(n-2,1)+h*(sigma*(w(n-2,2)-w(n-2,1)));
    w(n+1,2)=0*w(n,2)+1*w(n-2,2)+h*(w(n-2,1)*(rho-w(n-2,3))-w(n-2,2));
    w(n+1,3)=0*w(n,3)+1*w(n-2,3)+h*(w(n-2,1)*w(n-2,2)-beta*w(n-2,3));     
end
title('Lorenz system for explicit method')
plot3(w(:,1),w(:,2),w(:,3))
xlabel('x')
ylabel('y')
zlabel('z')
end