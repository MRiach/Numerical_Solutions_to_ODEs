%This function computes the global error for the Euler method
%used in question 1 
function globalerror = EulerGlobalError(tinit,tfinal,initcond,h)


 numsol = EulerMethod(tinit,tfinal,initcond,h);
 n = floor((tfinal-tinit)/h);
 %compute the global error by finding the difference between the official
 %value and the final value of the numerical solution
 globalerror = numsol(n+1,1) - ExactSolution(tfinal);
 globalerror = abs(globalerror);
 
end

