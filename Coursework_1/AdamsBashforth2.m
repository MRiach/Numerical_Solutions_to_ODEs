%This function implements the AB (2) Method using the f given in
%Exercise 3
function numsol = AdamsBashforth2(tinit,tfinal,initcond,h)

 n = (tfinal-tinit)/h;
 numsol = zeros(n+1,length(initcond));
 %Establish initial condition
 numsol(1,:) = initcond ;
 t = tinit:h:tfinal ;
 % Euler method is performed for the 2nd row since Adams-Bashforth(2) requires 2 rows of data before it can be
 % performed
 
 %Euler Method
 
 numsol(2,:) = numsol(1,:)+ h*[numsol(1,2),0.3*cos(t(1))-0.05*numsol(1,2)-numsol(1,1)^3];
 
 %Adams Bashforth 2 
 
 for i = 3:n+1
   numsol(i,:) = numsol(i-1,:) + 3*h/2*[numsol(i-1,2),0.3*cos(t(i-1))-0.05*numsol(i-1,2)-numsol(i-1,1)^3]...
       - h/2*[numsol(i-2,2),0.3*cos(t(i-2))-0.05*numsol(i-2,2)-numsol(i-2,1)^3];
 end 
 
end
