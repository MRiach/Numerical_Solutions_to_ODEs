%This function implements the Euler Method using the f given in Exercise 1
function numsol = EulerMethod(tinit,tfinal,initcond,h)
 
 n = floor((tfinal-tinit)/h);
 numsol = zeros(n+1,length(initcond));
 %Establish initial condition
 numsol(1,:) = initcond ;
 t = tinit:h:tfinal ;
 
 %Carry out Euler Method
 for i = 2:n+1
    numsol(i,:) = numsol(i-1,:)+ h*[numsol(i-1,2),t(i-1)^2-numsol(i-1,2)-4*numsol(i-1,1)];
 end 
 
end