%This function computes the global error for the AB3 Method
%used in question 5 
function globalerror = AB3GlobalError(tinit,tfinal,initcond,h)

 globalerror = zeros(1,2);
 numsol = AdamsBashforth3(tinit,tfinal,initcond,h);
 n = floor((tfinal-tinit)/h);
 %compute the global error by finding the difference between the official
 %value and the final value of the numerical solution
 globalerror(1) = numsol(n+1,1) - ExactSolutionx1(tfinal);
 globalerror(2) = numsol(n+1,2) - ExactSolutionx2(tfinal);
 globalerror = abs(globalerror);
 
end