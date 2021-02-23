%This function implements the AB (3) Method using the f given in
%Exercise 4
function numsol = AdamsBashforth3(tinit,tfinal,initcond,h)

 n = (tfinal-tinit)/h;
 numsol = zeros(n+1,length(initcond));
 %Establish initial condition
 numsol(1,:) = initcond ;
 % Euler method is performed for the 2nd row and AB(2) for the 3rd row of the numerical
 % solution, since Adams-Bashforth(3) requires 3 rows of data before it can be
 % performed
 
 %Euler Method
 
 numsol(2,:) = numsol(1,:)+ h*[10*(numsol(1,2)-0.1*numsol(1,1)^3+0.2*numsol(1,1)),numsol(1,1)-numsol(1,2)+numsol(1,3),-15*numsol(1,2)-0.01*numsol(1,3)];
 
 %Adams Bashforth 2 
 numsol(3,:) = numsol(2,:) + 3*h/2*[10*(numsol(2,2)-0.1*numsol(2,1)^3+0.2*numsol(2,1)),numsol(2,1)-numsol(2,2)+numsol(2,3),-15*numsol(2,2)-0.01*numsol(2,3)]...
     - h/2*[10*(numsol(1,2)-0.1*numsol(1,1)^3+0.2*numsol(1,1)),numsol(1,1)-numsol(1,2)+numsol(1,3),-15*numsol(1,2)-0.01*numsol(1,3)];
 
 %Adams Bashforth 3 
 for i = 4:n+1
    numsol(i,:) = numsol(i-1,:)+ 23*h/12*[10*(numsol(i-1,2)-0.1*numsol(i-1,1)^3+0.2*numsol(i-1,1)),numsol(i-1,1)-numsol(i-1,2)+numsol(i-1,3),-15*numsol(i-1,2)-0.01*numsol(i-1,3)]...
        - 4*h/3*[10*(numsol(i-2,2)-0.1*numsol(i-2,1)^3+0.2*numsol(i-2,1)),numsol(i-2,1)-numsol(i-2,2)+numsol(i-2,3),-15*numsol(i-2,2)-0.01*numsol(i-2,3)]...
        +5*h/12*[10*(numsol(i-3,2)-0.1*numsol(i-3,1)^3+0.2*numsol(i-3,1)),numsol(i-3,1)-numsol(i-3,2)+numsol(i-3,3),-15*numsol(i-3,2)-0.01*numsol(i-3,3)];
 end 
end


