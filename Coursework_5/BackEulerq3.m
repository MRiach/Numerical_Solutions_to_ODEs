 %Backward Euler applied to modified linear BVP in exercise 3

function numsol = BackEulerq3(T,L,dx,dt)

M = L/dx;
N = T/dt;

numsol = zeros(M+1,N+1);

X = 0:dx:L;
T = 0:dt:T;

%initial condition (all zeros)
numsol(:,1) = X*0;

%create matrix A for solving set of implicit equations and invert it 
%to produce solutions at next time step 

a = 1+dt*10+2*0.1*dt/(dx^2);
b = -0.1*dt/(dx^2);

A = diag(a*ones(M,1))+diag(b*ones(M-1,1),1)+diag(b*ones(M-1,1),-1);
A(1,M) = b;
A(M,1) = b;

for j = 2:N+1
    numsol(1:M,j) = inv(A)*(numsol(1:M,j-1)+dt*cos(2*pi*X(1:M)'/L)...
        *(1+0.1*T(j)*(2*pi/L)^2+10*T(j)));
end

%make final row of numsol the same as the first row so that boundary
%conditions are respected 
numsol(M+1,:) = numsol(1,:);

 
end