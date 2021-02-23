%Backward Euler applied to linear BVP in exercise 1 & 2 

function numsol = BackEulerq1(T,L,dx,dt)

M = L/dx;
N = T/dt;

numsol = zeros(M+1,N+1);

X0 = 0:dx:L;
X0 = exp(-L^2*(X0-L/2).^2);

%initial condition
numsol(:,1) = X0;

%create matrix A for solving set of implicit equations and invert it 
%to produce solutions at next time step 

a = 1+dt*10+2*0.1*dt/(dx^2);
b = -0.1*dt/(dx^2);

A = diag(a*ones(M,1))+diag(b*ones(M-1,1),1)+diag(b*ones(M-1,1),-1);
A(1,M) = b;
A(M,1) = b;

for j = 2:N+1
    numsol(1:M,j) = inv(A)*(numsol(1:M,j-1)+dt*10);
end

%make final row of numsol the same as the first row so that boundary
%conditions are respected 
numsol(M+1,:) = numsol(1,:);

 
end