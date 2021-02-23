%LMM method used in question 5b
function numsol = LMMq5b(tinit,tfinal,initcond,h)

 n = floor((tfinal-tinit)/h);
 numsol = zeros(n+1,length(initcond));
 %Establish initial condition
 numsol(1,:) = initcond;
 t = tinit:h:tfinal;
 % Euler method is performed for the 2nd row 
 
 %Euler Method
 
 numsol(2,:) = numsol(1,:)+ h*f1(numsol(1,1),numsol(1,2),numsol(1,3));
 
 %Adams Bashforth 2 and Adams Moulton 3 predictor/corrector  
 for i = 3:n+1
  %predictor   
  x_hat = numsol(i-1,:) + h/2*(3*f1(numsol(i-1,1),numsol(i-1,2),numsol(i-1,3))...
          -f1(numsol(i-2,1),numsol(i-2,2),numsol(i-2,3)));
  
  %corrector
  numsol(i,:) = numsol(i-1,:) + h*(5/12*f1(x_hat(1),x_hat(2),x_hat(3))...
                +2/3*f1(numsol(i-1,1),numsol(i-1,2),numsol(i-1,3))...
                 -1/12*f1(numsol(i-2,1),numsol(i-2,2),numsol(i-2,3)));
 end
 
end
