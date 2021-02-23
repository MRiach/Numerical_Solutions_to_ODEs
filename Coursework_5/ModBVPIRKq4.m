%IRK applied to modified nonlinear BVP in exercise 4

%m denotes the number of iterations used in the Newton method


function numsol = ModBVPIRKq4(T,L,dx,dt,v,m)

M = L/dx;
N = T/dt;

numsol = zeros(M+1,N+1);

X = 0:dx:L;
T = 0:dt:T;

%initial condition (all zeros)
numsol(:,1) = zeros(M+1,1);

%Input entries in Jacobian, J, which are constant for all iterations

J = zeros(M+1,M+1);
J(1,2) = -2*0.1*dt/(dx^2);
J(M+1,M) = -2*0.1*dt/(dx^2);
for i = 2:M
    J(i,i-1) = -0.1*dt/(dx^2)-v*dt/(2*dx);
    J(i,i+1) = -0.1*dt/(dx^2)+v*dt/(2*dx);
end 

%Implement 4 stage IRK 

for i = 2:N+1
    
   k1 = zeros(M+1,1);
   k2 = zeros(M+1,1);
   k3 = zeros(M+1,1);
   k4 = zeros(M+1,1);
   
   %k1
   
   u = numsol(:,i-1);
   k1(1) = 0.1/(dx^2)*(2*u(2)-2*u(1))+u(1)^2*(1-u(1))+sin(pi*X(1)/L+pi/2)+...
       0.1*pi^2/L^2*T(i-1)*sin(pi*X(1)/L+pi/2)+v*pi/L*T(i-1)*cos(pi*X(1)/L+pi/2)...
       -T(i-1)^2*(sin(pi*X(1)/L+pi/2))^2*(1-T(i-1)*sin(pi*X(1)/L+pi/2));
   k1(M+1) = 0.1/(dx^2)*(2*u(M)-2*u(M+1))+u(M+1)^2*(1-u(M+1))+sin(pi*X(M+1)/L+pi/2)+...
       0.1*pi^2/L^2*T(i-1)*sin(pi*X(M+1)/L+pi/2)+v*pi/L*T(i-1)*cos(pi*X(M+1)/L+pi/2)...
       -T(i-1)^2*(sin(pi*X(M+1)/L+pi/2))^2*(1-T(i-1)*sin(pi*X(M+1)/L+pi/2));
   for j = 2:M
       k1(j) = 0.1/(dx^2)*(u(j+1)-2*u(j)+u(j-1))-v/(2*dx)*(u(j+1)-u(j-1))+u(j)^2*(1-u(j))+sin(pi*X(j)/L+pi/2)+...
       0.1*pi^2/L^2*T(i-1)*sin(pi*X(j)/L+pi/2)+v*pi/L*T(i-1)*cos(pi*X(j)/L+pi/2)...
       -T(i-1)^2*(sin(pi*X(j)/L+pi/2))^2*(1-T(i-1)*sin(pi*X(j)/L+pi/2));      
   end   
   
   %k2
   
   u = numsol(:,i-1)+dt/3*k1;
   k2(1) = 0.1/(dx^2)*(2*u(2)-2*u(1))+u(1)^2*(1-u(1))+sin(pi*X(1)/L+pi/2)+...
       0.1*pi^2/L^2*(T(i-1)+dt/3)*sin(pi*X(1)/L+pi/2)+v*pi/L*(T(i-1)+dt/3)*cos(pi*X(1)/L+pi/2)...
       -(T(i-1)+dt/3)^2*(sin(pi*X(1)/L+pi/2))^2*(1-(T(i-1)+dt/3)*sin(pi*X(1)/L+pi/2));
   k2(M+1) = 0.1/(dx^2)*(2*u(M)-2*u(M+1))+u(M+1)^2*(1-u(M+1))+sin(pi*X(M+1)/L+pi/2)+...
       0.1*pi^2/L^2*(T(i-1)+dt/3)*sin(pi*X(M+1)/L+pi/2)+v*pi/L*(T(i-1)+dt/3)*cos(pi*X(M+1)/L+pi/2)...
       -(T(i-1)+dt/3)^2*(sin(pi*X(M+1)/L+pi/2))^2*(1-(T(i-1)+dt/3)*sin(pi*X(M+1)/L+pi/2));
   for j = 2:M
       k2(j) = 0.1/(dx^2)*(u(j+1)-2*u(j)+u(j-1))-v/(2*dx)*(u(j+1)-u(j-1))+u(j)^2*(1-u(j))+sin(pi*X(j)/L+pi/2)+...
       0.1*pi^2/L^2*(T(i-1)+dt/3)*sin(pi*X(j)/L+pi/2)+v*pi/L*(T(i-1)+dt/3)*cos(pi*X(j)/L+pi/2)...
       -(T(i-1)+dt/3)^2*(sin(pi*X(j)/L+pi/2))^2*(1-(T(i-1)+dt/3)*sin(pi*X(j)/L+pi/2));      
   end 
   
   %k3
   
   u = numsol(:,i-1)-dt/3*k1+dt*k2; 
   k3(1) = 0.1/(dx^2)*(2*u(2)-2*u(1))+u(1)^2*(1-u(1))+sin(pi*X(1)/L+pi/2)+...
       0.1*pi^2/L^2*(T(i-1)+2*dt/3)*sin(pi*X(1)/L+pi/2)+v*pi/L*(T(i-1)+2*dt/3)*cos(pi*X(1)/L+pi/2)...
       -(T(i-1)+2*dt/3)^2*(sin(pi*X(1)/L+pi/2))^2*(1-(T(i-1)+2*dt/3)*sin(pi*X(1)/L+pi/2));
   k3(M+1) = 0.1/(dx^2)*(2*u(M)-2*u(M+1))+u(M+1)^2*(1-u(M+1))+sin(pi*X(M+1)/L+pi/2)+...
       0.1*pi^2/L^2*(T(i-1)+2*dt/3)*sin(pi*X(M+1)/L+pi/2)+v*pi/L*(T(i-1)+2*dt/3)*cos(pi*X(M+1)/L+pi/2)...
       -(T(i-1)+2*dt/3)^2*(sin(pi*X(M+1)/L+pi/2))^2*(1-(T(i-1)+2*dt/3)*sin(pi*X(M+1)/L+pi/2));
   for j = 2:M
       k3(j) = 0.1/(dx^2)*(u(j+1)-2*u(j)+u(j-1))-v/(2*dx)*(u(j+1)-u(j-1))+u(j)^2*(1-u(j))+sin(pi*X(j)/L+pi/2)+...
       0.1*pi^2/L^2*(T(i-1)+2*dt/3)*sin(pi*X(j)/L+pi/2)+v*pi/L*(T(i-1)+2*dt/3)*cos(pi*X(j)/L+pi/2)...
       -(T(i-1)+2*dt/3)^2*(sin(pi*X(j)/L+pi/2))^2*(1-(T(i-1)+2*dt/3)*sin(pi*X(j)/L+pi/2));      
   end 
   %initialise guess for k4
   
   u = numsol(:,i-1)+2*dt*k1-dt*k3;
   k4(1) = 0.1/(dx^2)*(2*u(2)-2*u(1))+u(1)^2*(1-u(1))+sin(pi*X(1)/L+pi/2)+...
       0.1*pi^2/L^2*(T(i-1)+dt)*sin(pi*X(1)/L+pi/2)+v*pi/L*(T(i-1)+dt)*cos(pi*X(1)/L+pi/2)...
       -(T(i-1)+dt)^2*(sin(pi*X(1)/L+pi/2))^2*(1-(T(i-1)+dt)*sin(pi*X(1)/L+pi/2));
   k4(M+1) = 0.1/(dx^2)*(2*u(M)-2*u(M+1))+u(M+1)^2*(1-u(M+1))+sin(pi*X(M+1)/L+pi/2)+...
       0.1*pi^2/L^2*(T(i-1)+dt)*sin(pi*X(M+1)/L+pi/2)+v*pi/L*(T(i-1)+dt)*cos(pi*X(M+1)/L+pi/2)...
       -(T(i-1)+dt)^2*(sin(pi*X(M+1)/L+pi/2))^2*(1-(T(i-1)+dt)*sin(pi*X(M+1)/L+pi/2));
   for j = 2:M
       k4(j) = 0.1/(dx^2)*(u(j+1)-2*u(j)+u(j-1))-v/(2*dx)*(u(j+1)-u(j-1))+u(j)^2*(1-u(j))+sin(pi*X(j)/L+pi/2)+...
       0.1*pi^2/L^2*(T(i-1)+dt)*sin(pi*X(j)/L+pi/2)+v*pi/L*(T(i-1)+dt)*cos(pi*X(j)/L+pi/2)...
       -(T(i-1)+dt)^2*(sin(pi*X(j)/L+pi/2))^2*(1-(T(i-1)+dt)*sin(pi*X(j)/L+pi/2));      
   end 
   
   %configure Jacobian and value of F at 0th iteration 
   
   %Jacobian
   
   a = numsol(:,i-1)+dt*k1-dt*k3;
   for j = 1:M+1
       J(j,j) = 1+2*0.1*dt/(dx^2)+(2*dt^2*k4(j)+2*dt*a(j))*(1-a(j)-dt*k4(j))-...
           dt*(a(j)^2+2*a(j)*dt*k4(j)+dt^2*k4(j)^2);
   end
   
   %F
   
   F = zeros(M+1,1);
   u = numsol(:,i-1)+dt*k1-dt*k3+dt*k4;
   F(1) = 0.1/(dx^2)*(2*u(2)-2*u(1))+u(1)^2*(1-u(1))+sin(pi*X(1)/L+pi/2)+...
       0.1*pi^2/L^2*(T(i-1)+dt)*sin(pi*X(1)/L+pi/2)+v*pi/L*(T(i-1)+dt)*cos(pi*X(1)/L+pi/2)...
       -(T(i-1)+dt)^2*(sin(pi*X(1)/L+pi/2))^2*(1-(T(i-1)+dt)*sin(pi*X(1)/L+pi/2));
   F(M+1) = 0.1/(dx^2)*(2*u(M)-2*u(M+1))+u(M+1)^2*(1-u(M+1))+sin(pi*X(M+1)/L+pi/2)+...
       0.1*pi^2/L^2*(T(i-1)+dt)*sin(pi*X(M+1)/L+pi/2)+v*pi/L*(T(i-1)+dt)*cos(pi*X(M+1)/L+pi/2)...
       -(T(i-1)+dt)^2*(sin(pi*X(M+1)/L+pi/2))^2*(1-(T(i-1)+dt)*sin(pi*X(M+1)/L+pi/2));
   for j = 2:M
       F(j) = 0.1/(dx^2)*(u(j+1)-2*u(j)+u(j-1))-v/(2*dx)*(u(j+1)-u(j-1))+u(j)^2*(1-u(j))+sin(pi*X(j)/L+pi/2)+...
       0.1*pi^2/L^2*(T(i-1)+dt)*sin(pi*X(j)/L+pi/2)+v*pi/L*(T(i-1)+dt)*cos(pi*X(j)/L+pi/2)...
       -(T(i-1)+dt)^2*(sin(pi*X(j)/L+pi/2))^2*(1-(T(i-1)+dt)*sin(pi*X(j)/L+pi/2));      
   end 
   F = k4-F;
   
   
   for k = 1:m
       
     k4a = k4 - inv(J)*F;   
     
     %implement stopping condition
     
     if norm(k4a-k4)/norm(k4a) < 1e-5
        k4 = k4a; 
        break
     end
     
     %move onto next iteration, updating Jacobian and F prior to proceeding
     %to next step
     
     %update k4 as new iteration
     
     k4 = k4a;
     
     %Update entries of Jacobian (only diagonal elements are affected)
     
     for j = 1:M+1
       J(j,j) = 1+2*0.1*dt/(dx^2)+(2*dt^2*k4(j)+2*dt*a(j))*(1-a(j)-dt*k4(j))-...
           dt*(a(j)^2+2*a(j)*dt*k4(j)+dt^2*k4(j)^2);
     end
     
     %Evaluate F at new iteration 
     F = zeros(M+1,1);
     u = numsol(:,i-1)+dt*k1-dt*k3+dt*k4;
     F(1) = 0.1/(dx^2)*(2*u(2)-2*u(1))+u(1)^2*(1-u(1))+sin(pi*X(1)/L+pi/2)+...
       0.1*pi^2/L^2*(T(i-1)+dt)*sin(pi*X(1)/L+pi/2)+v*pi/L*(T(i-1)+dt)*cos(pi*X(1)/L+pi/2)...
       -(T(i-1)+dt)^2*(sin(pi*X(1)/L+pi/2))^2*(1-(T(i-1)+dt)*sin(pi*X(1)/L+pi/2));
     F(M+1) = 0.1/(dx^2)*(2*u(M)-2*u(M+1))+u(M+1)^2*(1-u(M+1))+sin(pi*X(M+1)/L+pi/2)+...
       0.1*pi^2/L^2*(T(i-1)+dt)*sin(pi*X(M+1)/L+pi/2)+v*pi/L*(T(i-1)+dt)*cos(pi*X(M+1)/L+pi/2)...
       -(T(i-1)+dt)^2*(sin(pi*X(M+1)/L+pi/2))^2*(1-(T(i-1)+dt)*sin(pi*X(M+1)/L+pi/2));
     for j = 2:M
       F(j) = 0.1/(dx^2)*(u(j+1)-2*u(j)+u(j-1))-v/(2*dx)*(u(j+1)-u(j-1))+u(j)^2*(1-u(j))+sin(pi*X(j)/L+pi/2)+...
       0.1*pi^2/L^2*(T(i-1)+dt)*sin(pi*X(j)/L+pi/2)+v*pi/L*(T(i-1)+dt)*cos(pi*X(j)/L+pi/2)...
       -(T(i-1)+dt)^2*(sin(pi*X(j)/L+pi/2))^2*(1-(T(i-1)+dt)*sin(pi*X(j)/L+pi/2));      
     end 
     F = k4-F;
     
     %if stop condition is never reached, assign solution as value from
     %last iteration     
     
     if j == m
         k4 = k4a; 
     end
     
   end  
   
   numsol(:,i) = numsol(:,i-1)+dt/8*(k1+3*k2+3*k3+k4);
end

 
end