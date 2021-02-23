%This function computes the global error for the TS(3) method
%used in question 2 
function globalerror = TaylorSeriesGlobalError(tinit,tfinal,initcond,h)


 numsol = TaylorSeriesMethod(tinit,tfinal,initcond,h);
 n = (tfinal-tinit)/h;
 %compute the global error by finding the difference between the official
 %value and the final value of the numerical solution
 globalerror = numsol(n+1,1) - ExactSolution1(tfinal);
 globalerror = abs(globalerror);
 
end

