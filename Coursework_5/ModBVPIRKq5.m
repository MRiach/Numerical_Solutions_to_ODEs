%IRK applied to modified nonlinear BVP in exercise 5

%m denotes the number of iterations used in the fixed point iterations


function numsol = ModBVPIRKq5(T,L,dx,dt,m)

M = L/dx;
N = T/dt;

numsol = zeros(M+1,N+1);

X = 0:dx:L;
T = 0:dt:T;

%initial condition 

numsol(:,1) = zeros(M+1,1);

%initialise matrix A and compute its inverse for finding u_t which is 
%constant throughout the algorithm

a = 1+2/(dx^2);
b = -1/(dx^2);

A = diag(a*ones(M+1,1))+diag(b*ones(M,1),1)+diag(b*ones(M,1),-1);
A(1,M) = b;
A(M+1,2) = b;
Ainv = inv(A);

%Implement 4 stage IRK 

for i = 2:N+1
    
   k1 = zeros(M+1,1);
   k2 = zeros(M+1,1);
   k3 = zeros(M+1,1);
   k4 = zeros(M+1,1);
   
   %k1
   
   %1st component and M+1 th component are equivalent by boundary
   %conditions
   
   u = numsol(:,i-1);
   
   k1(1) = -(u(1)-(u(2)-2*u(1)+u(M))/(dx)^2)*(u(2)-u(M))/(2*dx)+...
       u(1)*(u(2)-2*u(1)+u(M))/(dx)^3-u(M)*(u(1)-2*u(M)+u(M-1))/(dx)^3 ... 
        -(u(1)^2-u(M)^2)/dx+(1+(2*pi/L)^2)*cos(2*pi*X(1)/L)-3*pi/L*(1+(2*pi/L)^2)*T(i-1)^2*sin(4*pi*X(1)/L);
    
   k1(2) = -(u(2)-(u(3)-2*u(2)+u(1))/(dx)^2)*(u(3)-u(1))/(2*dx)+...
       u(2)*(u(3)-2*u(2)+u(1))/(dx)^3-u(1)*(u(2)-2*u(1)+u(M))/(dx)^3 ... 
        -(u(2)^2-u(1)^2)/dx+(1+(2*pi/L)^2)*cos(2*pi*X(2)/L)-3*pi/L*(1+(2*pi/L)^2)*T(i-1)^2*sin(4*pi*X(2)/L);
   
   k1(M+1) = -(u(1)-(u(2)-2*u(1)+u(M))/(dx)^2)*(u(2)-u(M))/(2*dx)+...
       u(1)*(u(2)-2*u(1)+u(M))/(dx)^3-u(M)*(u(1)-2*u(M)+u(M-1))/(dx)^3 ... 
        -(u(1)^2-u(M)^2)/dx+(1+(2*pi/L)^2)*cos(2*pi*X(M+1)/L)-3*pi/L*(1+(2*pi/L)^2)*T(i-1)^2*sin(4*pi*X(M+1)/L);
   
   
   for j = 3:M 
       k1(j) = -(u(j)-(u(j+1)-2*u(j)+u(j-1))/(dx)^2)*(u(j+1)-u(j-1))/(2*dx)+...
       u(j)*(u(j+1)-2*u(j)+u(j-1))/(dx)^3-u(j-1)*(u(j)-2*u(j-1)+u(j-2))/(dx)^3 ... 
        -(u(j)^2-u(j-1)^2)/dx+(1+(2*pi/L)^2)*cos(2*pi*X(j)/L)-3*pi/L*(1+(2*pi/L)^2)*T(i-1)^2*sin(4*pi*X(j)/L);      
   end    
   
   %Invert to isolate u_t on LHS
   
   k1 = Ainv*k1;
   
   %k2
   
   u = numsol(:,i-1)+dt/3*k1;
   
   k2(1) = -(u(1)-(u(2)-2*u(1)+u(M))/(dx)^2)*(u(2)-u(M))/(2*dx)+...
       u(1)*(u(2)-2*u(1)+u(M))/(dx)^3-u(M)*(u(1)-2*u(M)+u(M-1))/(dx)^3 ... 
        -(u(1)^2-u(M)^2)/dx+(1+(2*pi/L)^2)*cos(2*pi*X(1)/L)-3*pi/L*(1+(2*pi/L)^2)*(T(i-1)+dt/3)^2*sin(4*pi*X(1)/L);
    
   k2(2) = -(u(2)-(u(3)-2*u(2)+u(1))/(dx)^2)*(u(3)-u(1))/(2*dx)+...
       u(2)*(u(3)-2*u(2)+u(1))/(dx)^3-u(1)*(u(2)-2*u(1)+u(M))/(dx)^3 ... 
        -(u(2)^2-u(1)^2)/dx+(1+(2*pi/L)^2)*cos(2*pi*X(2)/L)-3*pi/L*(1+(2*pi/L)^2)*(T(i-1)+dt/3)^2*sin(4*pi*X(2)/L);
   
   k2(M+1) = -(u(1)-(u(2)-2*u(1)+u(M))/(dx)^2)*(u(2)-u(M))/(2*dx)+...
       u(1)*(u(2)-2*u(1)+u(M))/(dx)^3-u(M)*(u(1)-2*u(M)+u(M-1))/(dx)^3 ... 
        -(u(1)^2-u(M)^2)/dx+(1+(2*pi/L)^2)*cos(2*pi*X(M+1)/L)-3*pi/L*(1+(2*pi/L)^2)*(T(i-1)+dt/3)^2*sin(4*pi*X(M+1)/L);
   
   
   for j = 3:M
       k2(j) = -(u(j)-(u(j+1)-2*u(j)+u(j-1))/(dx)^2)*(u(j+1)-u(j-1))/(2*dx)+...
       u(j)*(u(j+1)-2*u(j)+u(j-1))/(dx)^3-u(j-1)*(u(j)-2*u(j-1)+u(j-2))/(dx)^3 ... 
        -(u(j)^2-u(j-1)^2)/dx+(1+(2*pi/L)^2)*cos(2*pi*X(j)/L)-3*pi/L*(1+(2*pi/L)^2)*(T(i-1)+dt/3)^2*sin(4*pi*X(j)/L);      
   end    
   
   k2 = Ainv*k2;
      
   %k3
   
   u = numsol(:,i-1)-dt/3*k1+dt*k2;
   
   
   k3(1) = -(u(1)-(u(2)-2*u(1)+u(M))/(dx)^2)*(u(2)-u(M))/(2*dx)+...
       u(1)*(u(2)-2*u(1)+u(M))/(dx)^3-u(M)*(u(1)-2*u(M)+u(M-1))/(dx)^3 ... 
        -(u(1)^2-u(M)^2)/dx+(1+(2*pi/L)^2)*cos(2*pi*X(1)/L)-3*pi/L*(1+(2*pi/L)^2)*(T(i-1)+2*dt/3)^2*sin(4*pi*X(1)/L);
    
   k3(2) = -(u(2)-(u(3)-2*u(2)+u(1))/(dx)^2)*(u(3)-u(1))/(2*dx)+...
       u(2)*(u(3)-2*u(2)+u(1))/(dx)^3-u(1)*(u(2)-2*u(1)+u(M))/(dx)^3 ... 
        -(u(2)^2-u(1)^2)/dx+(1+(2*pi/L)^2)*cos(2*pi*X(2)/L)-3*pi/L*(1+(2*pi/L)^2)*(T(i-1)+2*dt/3)^2*sin(4*pi*X(2)/L);
   
   k3(M+1) = -(u(1)-(u(2)-2*u(1)+u(M))/(dx)^2)*(u(2)-u(M))/(2*dx)+...
       u(1)*(u(2)-2*u(1)+u(M))/(dx)^3-u(M)*(u(1)-2*u(M)+u(M-1))/(dx)^3 ... 
        -(u(1)^2-u(M)^2)/dx+(1+(2*pi/L)^2)*cos(2*pi*X(M+1)/L)-3*pi/L*(1+(2*pi/L)^2)*(T(i-1)+2*dt/3)^2*sin(4*pi*X(M+1)/L);
   
   
   for j = 3:M
       k3(j) = -(u(j)-(u(j+1)-2*u(j)+u(j-1))/(dx)^2)*(u(j+1)-u(j-1))/(2*dx)+...
       u(j)*(u(j+1)-2*u(j)+u(j-1))/(dx)^3-u(j-1)*(u(j)-2*u(j-1)+u(j-2))/(dx)^3 ... 
        -(u(j)^2-u(j-1)^2)/dx+(1+(2*pi/L)^2)*cos(2*pi*X(j)/L)-3*pi/L*(1+(2*pi/L)^2)*(T(i-1)+2*dt/3)^2*sin(4*pi*X(j)/L);      
   end    
   
   k3 = Ainv*k3;

   %initialise guess for k4
   
   u = numsol(:,i-1)+2*dt*k1-dt*k3;
   
   k4(1) = -(u(1)-(u(2)-2*u(1)+u(M))/(dx)^2)*(u(2)-u(M))/(2*dx)+...
       u(1)*(u(2)-2*u(1)+u(M))/(dx)^3-u(M)*(u(1)-2*u(M)+u(M-1))/(dx)^3 ... 
        -(u(1)^2-u(M)^2)/dx+(1+(2*pi/L)^2)*cos(2*pi*X(1)/L)-3*pi/L*(1+(2*pi/L)^2)*(T(i-1)+dt)^2*sin(4*pi*X(1)/L);
    
   k4(2) = -(u(2)-(u(3)-2*u(2)+u(1))/(dx)^2)*(u(3)-u(1))/(2*dx)+...
       u(2)*(u(3)-2*u(2)+u(1))/(dx)^3-u(1)*(u(2)-2*u(1)+u(M))/(dx)^3 ... 
        -(u(2)^2-u(1)^2)/dx+(1+(2*pi/L)^2)*cos(2*pi*X(2)/L)-3*pi/L*(1+(2*pi/L)^2)*(T(i-1)+dt)^2*sin(4*pi*X(2)/L);
   
   k4(M+1) = -(u(1)-(u(2)-2*u(1)+u(M))/(dx)^2)*(u(2)-u(M))/(2*dx)+...
       u(1)*(u(2)-2*u(1)+u(M))/(dx)^3-u(M)*(u(1)-2*u(M)+u(M-1))/(dx)^3 ... 
        -(u(1)^2-u(M)^2)/dx+(1+(2*pi/L)^2)*cos(2*pi*X(M+1)/L)-3*pi/L*(1+(2*pi/L)^2)*(T(i-1)+dt)^2*sin(4*pi*X(M+1)/L);
   
   
   for j = 3:M
       k4(j) = -(u(j)-(u(j+1)-2*u(j)+u(j-1))/(dx)^2)*(u(j+1)-u(j-1))/(2*dx)+...
       u(j)*(u(j+1)-2*u(j)+u(j-1))/(dx)^3-u(j-1)*(u(j)-2*u(j-1)+u(j-2))/(dx)^3 ... 
        -(u(j)^2-u(j-1)^2)/dx+(1+(2*pi/L)^2)*cos(2*pi*X(j)/L)-3*pi/L*(1+(2*pi/L)^2)*(T(i-1)+dt)^2*sin(4*pi*X(j)/L);      
   end    
   
   k4 = Ainv*k4;
   
   
   for k = 1:m
     
     k4a = zeros(M+1,1);
     u = numsol(:,i-1)+dt*k1-dt*k3+dt*k4;  
        
     k4a(1) = -(u(1)-(u(2)-2*u(1)+u(M))/(dx)^2)*(u(2)-u(M))/(2*dx)+...
       u(1)*(u(2)-2*u(1)+u(M))/(dx)^3-u(M)*(u(1)-2*u(M)+u(M-1))/(dx)^3 ... 
        -(u(1)^2-u(M)^2)/dx+(1+(2*pi/L)^2)*cos(2*pi*X(1)/L)-3*pi/L*(1+(2*pi/L)^2)*(T(i-1)+dt)^2*sin(4*pi*X(1)/L);
    
     k4a(2) = -(u(2)-(u(3)-2*u(2)+u(1))/(dx)^2)*(u(3)-u(1))/(2*dx)+...
       u(2)*(u(3)-2*u(2)+u(1))/(dx)^3-u(1)*(u(2)-2*u(1)+u(M))/(dx)^3 ... 
        -(u(2)^2-u(1)^2)/dx+(1+(2*pi/L)^2)*cos(2*pi*X(2)/L)-3*pi/L*(1+(2*pi/L)^2)*(T(i-1)+dt)^2*sin(4*pi*X(2)/L);
   
     k4a(M+1) = -(u(1)-(u(2)-2*u(1)+u(M))/(dx)^2)*(u(2)-u(M))/(2*dx)+...
       u(1)*(u(2)-2*u(1)+u(M))/(dx)^3-u(M)*(u(1)-2*u(M)+u(M-1))/(dx)^3 ... 
        -(u(1)^2-u(M)^2)/dx+(1+(2*pi/L)^2)*cos(2*pi*X(M+1)/L)-3*pi/L*(1+(2*pi/L)^2)*(T(i-1)+dt)^2*sin(4*pi*X(M+1)/L);
   
   
     for j = 3:M
       k4a(j) = -(u(j)-(u(j+1)-2*u(j)+u(j-1))/(dx)^2)*(u(j+1)-u(j-1))/(2*dx)+...
       u(j)*(u(j+1)-2*u(j)+u(j-1))/(dx)^3-u(j-1)*(u(j)-2*u(j-1)+u(j-2))/(dx)^3 ... 
        -(u(j)^2-u(j-1)^2)/dx+(1+(2*pi/L)^2)*cos(2*pi*X(j)/L)-3*pi/L*(1+(2*pi/L)^2)*(T(i-1)+dt)^2*sin(4*pi*X(j)/L);      
     end  
     
     k4a = Ainv*k4a;
     %implement stopping condition
     
     if norm(k4a-k4)/norm(k4a) < 1e-16
        k4 = k4a; 
        break
     end
          
     %update k4 as new iteration
     
     k4 = k4a;
          
     %if stop condition is never reached, assign solution as value from
     %last iteration     
     
     if j == m
         k4 = k4a; 
     end
     
   end  
   
   numsol(:,i) = numsol(:,i-1)+dt/8*(k1+3*k2+3*k3+k4);
end

 
end