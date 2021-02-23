%LMM method used in question 1
function numsol = LMMq1(tinit,tfinal,initcond,h)

 n = floor((tfinal-tinit)/h);
 numsol = zeros(n+1,length(initcond));
 %Establish initial condition
 numsol(1,:) = initcond;
 t = tinit:h:tfinal;
 % Euler method is performed for the 2nd row and AB(2) for the 3rd row of the numerical
 % solution, since LMM(3) requires 3 rows of data before it can be
 % performed
 
 %Euler Method
 
 numsol(2,:) = numsol(1,:)+ h*f(numsol(1,1),numsol(1,2),t(1));
 
 %Adams Bashforth 2 
 numsol(3,:) = numsol(2,:) + h/2*(3*f(numsol(2,1),numsol(2,2),t(2))-f(numsol(1,1),numsol(1,2),t(1)));
 
 %LMM 3 step
 for i = 4:n+1
  numsol(i,:) = numsol(i-3,:)+3*h*f(numsol(i-3,1),numsol(i-3,2),t(i-3));
 end
 
end

