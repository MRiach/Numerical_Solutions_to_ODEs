%This function gives the exact solution for Exercise 1
function x = ExactSolution(t)

 x = exp(-0.5*t)*((75/(32*sqrt(15)))*sin(sqrt(15)*t/2)+(3/32)*cos(sqrt(15)*t/2))+t^2/4-t/8-3/32;
 
end

