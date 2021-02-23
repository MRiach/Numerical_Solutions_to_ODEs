%AB5 method used in question 1
function numsol = AB5(tinit,initcond,h)

 tfinal = h*10000;
 n = floor((tfinal-tinit)/h);
 numsol = zeros(n+1,length(initcond));
 %Establish initial condition
 numsol(1,:) = initcond;
 t = tinit:h:tfinal;
 %Euler Method, AB1, AB2, AB3 and AB4 are all carried out prior to AB5
 
 %Euler Method
 
 numsol(2,:) = numsol(1,:)+ h*f(numsol(1,1),numsol(1,2));
 
 %Adams Bashforth 2 
 numsol(3,:) = numsol(2,:) + h/2*(3*f(numsol(2,1),numsol(2,2))-f(numsol(1,1),numsol(1,2)));
 
 %Adams Bashforth 3 
 numsol(4,:) = numsol(3,:) + h/12*(23*f(numsol(3,1),numsol(3,2))-...
                    16*f(numsol(2,1),numsol(2,2))+5*f(numsol(1,1),numsol(1,2)));
 
 %Adams Bashforth 4 
 numsol(5,:) = numsol(4,:) + h/24*(55*f(numsol(4,1),numsol(4,2))-59*f(numsol(3,1),numsol(3,2))...
     +37*f(numsol(2,1),numsol(2,2))-9*f(numsol(1,1),numsol(1,2)));
 
 %Adams Bashforth 5 
 for i = 6:n+1
  numsol(i,:) = numsol(i-1,:)+h/720*(1901*f(numsol(i-1,1),numsol(i-1,2))-2774*f(numsol(i-2,1),numsol(i-2,2))...
     +2616*f(numsol(i-3,1),numsol(i-3,2))-1274*f(numsol(i-4,1),numsol(i-4,2))+251*f(numsol(i-5,1),numsol(i-5,2)));
 end
 
end