%Exact solution to q4

function sol = ExactSolutionq4(T,L,dx,dt)

M = L/dx;
N = T/dt;

sol = zeros(M+1,N+1);
X = 0:dx:L;
T = 0:dt:T;

for j = 1:N+1
    sol(:,j) = T(j)*sin(pi*X/L+pi/2);
end 


end