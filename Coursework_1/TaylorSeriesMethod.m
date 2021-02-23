%This function implements the Taylor Series 2 Method using the f given in
%Exercise 2
function numsol = TaylorSeriesMethod(tinit,tfinal,initcond,h)
 
 n = (tfinal-tinit)/h;
 numsol = zeros(n+1,length(initcond));
 %Establish initial condition
 numsol(1,:) = initcond ;
 t = tinit:h:tfinal ;
 
 %Carry out TS(3) Method
  for i = 2:n+1
    numsol(i,:) = numsol(i-1,:)+ h*[numsol(i-1,2),t(i-1)^2-3*numsol(i-1,2)-2*numsol(i-1,1)]...
        +h^2/2*[t(i-1)^2-3*numsol(i-1,2)-2*numsol(i-1,1),2*t(i-1)-3*t(i-1)^2+7*numsol(i-1,2)+6*numsol(i-1,1)]...
        +h^3/6*[2*t(i-1)-3*t(i-1)^2+7*numsol(i-1,2)+6*numsol(i-1,1),7*t(i-1)^2-6*t(i-1)+2-15*numsol(i-1,2)-14*numsol(i-1,1)];
 end 