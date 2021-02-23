%LMM method used in question 6
function numsol = LMMq6(tinit,tfinal,initcond,h,M)

 n = floor((tfinal-tinit)/h);
 numsol = zeros(n+1,length(initcond));
 alpha2 = 9/10;
 alpha1 = 1/10;
 beta3 = -151/240;
 beta2 = 923/240;
 beta1 = -757/240;
 beta0 = 83/80;
 
 %Establish initial condition
 numsol(1,:) = initcond;
 t = tinit:h:tfinal;
 % Euler method is performed for the 2nd row and AB(2) for the 3rd row of the numerical
 % solution, since implicit method requires 3 rows of data before it can be
 % performed
 
 %Euler Method
 
 numsol(2,:) = numsol(1,:)+ h*f1(numsol(1,1),numsol(1,2),numsol(1,3));
 
 %Adams Bashforth 2 
 numsol(3,:) = numsol(2,:) + h/2*(3*f1(numsol(2,1),numsol(2,2),numsol(2,3))-f1(numsol(1,1),numsol(1,2),numsol(1,3)));
 
 %Newton method to implement implicit LMM  
 for i = 4:n+1
  %Let initial guess be the value at the prior time step due to the continuity
  %of the system
  xj =  numsol(i-1,:);  
  
  %apply M steps of Newton method using Fprime which represents the
  %Jacobian
  for j = 1:M
     xj1 = xj - (inv(Fprime(xj(1),xj(2),xj(3),h,beta3))*(xj-(alpha2*numsol(i-1,:)...
         +alpha1*numsol(i-2,:)+h*(beta3*f1(xj(1),xj(2),xj(3))+beta2*f1(numsol(i-1,1),numsol(i-1,2),numsol(i-1,3))...
         +beta1*f1(numsol(i-2,1),numsol(i-2,2),numsol(i-2,3))+beta0*f1(numsol(i-3,1),numsol(i-3,2),numsol(i-3,3)))))')';
     %implement stopping condition
     if norm(xj1-xj)/norm(xj1) < 1e-5
        numsol(i,:) = xj1; 
        break
     end
     %move onto next iteration
     xj = xj1;
     %if stop condition is never reached, assign solution as value from
     %last iteration     
     if j ==M
         numsol(i,:) = xj1; 
     end
  end   
 end 
end