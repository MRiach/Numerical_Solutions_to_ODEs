%Implicit Runge Kutta method used in question 5 with projection
function numsol = IRKq5proj(tinit,tfinal,initcond,h,M)

 n = floor((tfinal-tinit)/h);
 numsol = zeros(n+1,length(initcond));
 %Establish initial condition
 numsol(1,:) = initcond;
 t = tinit:h:tfinal;
 
 %IRK 4 stage and projection
 for i = 2:n+1
   k1 = f5(numsol(i-1,1),numsol(i-1,2),numsol(i-1,3));
   k2 = f5(numsol(i-1,1)+h/3*k1(1),numsol(i-1,2)+h/3*k1(2),numsol(i-1,3)+h/3*k1(3));
   k3 = f5(numsol(i-1,1)-h/3*k1(1)+h*k2(1),numsol(i-1,2)-h/3*k1(2)+h*k2(2),numsol(i-1,3)-h/3*k1(3)+h*k2(3));
   %initialise guess for k4
   k4 = f5(numsol(i-1,1)+2*h*k1(1)-h*k3(1),numsol(i-1,2)+2*h*k1(2)-h*k3(2),numsol(i-1,3)+2*h*k1(3)-h*k3(3));
   a = numsol(i-1,:)+h*k1-h*k3;
   for j = 1:M
     k4a = k4 - (inv(Fprimeq5(k4(1),k4(2),k4(3),a(1),a(2),a(3),h))*(k4-f5(a(1)+h*k4(1),...
              a(2)+h*k4(2),a(3)+h*k4(3)))')';    
     %implement stopping condition
     if norm(k4a-k4)/norm(k4a) < 1e-5
        k4 = k4a; 
        break
     end
     %move onto next iteration
     k4 = k4a;
     %if stop condition is never reached, assign solution as value from
     %last iteration     
     if j ==M
         k4 = k4a; 
     end
   end  
   %project solution of next step according to IRK onto sphere 
   y = numsol(i-1,:)+h/8*(k1+3*k2+3*k3+k4);
   numsol(i,:) = [-1,3,0.5]+3*(y-[-1,3,0.5])/norm(y-[-1,3,0.5]);
 end
 
 
end
