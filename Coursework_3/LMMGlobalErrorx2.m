%This function computes the global error for the LMM
%used in question 1 to the IVP in question 2 - 2nd component
function globalerrorx2 = LMMGlobalErrorx2(tinit,tfinal,initcond,h)

 numsol = LMMq1(tinit,tfinal,initcond,h);
 t = tinit:h:tfinal;
 n = floor((tfinal-tinit)/h);
 %compute the global error by finding the difference between the official
 %value and the final value of the numerical solution
 
 ExactSolution = ExactSolutionx2(t);
 ExactSolution = reshape(ExactSolution,[n+1,1]);
 globalerrorx2 = numsol(:,2) - ExactSolution;
 globalerrorx2 = abs(globalerrorx2);
 
end