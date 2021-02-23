%Exact solution to q3

function sol = ExactSolutionq3(T,L,dx,dt)

M = L/dx;
N = T/dt;

sol = zeros(M+1,N+1);
X = 0:dx:L;
T = 0:dt:T;

for j = 1:N+1
    sol(:,j) = T(j)*cos(2*pi*X/L);
end 


end