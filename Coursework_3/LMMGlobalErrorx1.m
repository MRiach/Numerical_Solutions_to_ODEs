%This function computes the global error for the LMM
%used in question 1 to the IVP in question 2 - 1st component
function globalerrorx1 = LMMGlobalErrorx1(tinit,tfinal,initcond,h)

 numsol = LMMq1(tinit,tfinal,initcond,h);
 t = tinit:h:tfinal;
 n = floor((tfinal-tinit)/h);
 %compute the global error by finding the difference between the official
 %value and the final value of the numerical solution
 
 ExactSolution = ExactSolutionx1(t);
 ExactSolution = reshape(ExactSolution,[n+1,1]);
 globalerrorx1 = numsol(:,1) - ExactSolution;
 globalerrorx1 = abs(globalerrorx1);
 
end

%Plots of global errors for different h vs time
%X = LMMGlobalErrorx1(0,20,[2,3],0.0005); h = 1/2000
%X = LMMGlobalErrorx1(0,20,[2,3],1/1499); h = 1/1499

%Plot numerical solution and official solution too 