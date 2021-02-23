%Explicit Runge Kutta method used in question 4
function numsol = ERK(tinit,tfinal,initcond,h)

 n = floor((tfinal-tinit)/h);
 numsol = zeros(n+1,length(initcond));
 %Establish initial condition
 numsol(1,:) = initcond;
 t = tinit:h:tfinal;
 
 %ERK 4 stage
 for i = 2:n+1
   k1 = f4(numsol(i-1,1),numsol(i-1,2),numsol(i-1,3));
   k2 = f4(numsol(i-1,1)+h/2*k1(1),numsol(i-1,2)+h/2*k1(2),numsol(i-1,3)+h/2*k1(3));
   k3 = f4(numsol(i-1,1)+h/3*k1(1)+h/6*k2(1),numsol(i-1,2)+h/3*k1(2)+h/6*k2(2),numsol(i-1,3)+h/3*k1(3)+h/6*k2(3));
   k4 = f4(numsol(i-1,1)-2*h*k2(1)+3*h*k3(1),numsol(i-1,2)-2*h*k2(2)+3*h*k3(2),numsol(i-1,3)-2*h*k2(3)+3*h*k3(3));
   numsol(i,:) = numsol(i-1,:)+h/6*(k1-2*k2+6*k3+k4);
 end
 
end
