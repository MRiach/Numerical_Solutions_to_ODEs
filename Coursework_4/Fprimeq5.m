%This produces the Jacobian used in Newton method in q5
function matrix = Fprimeq5(x,y,z,a1,a2,a3,h)
 matrix = [1-2*h*(a2+h*y)*(a1+h*x)-0.87*h,-h*(a3+h*z-1+(a1+h*x)^2)...
     ,-h*(a2+h*y);h*(2*(a1+h*x)^2-(3*a3+3*h*z+1-(a1+h*x)^2)),1-0.87*h...
     ,-3*h*(a1+h*x);2*(a3+h*z)*(h*a2+y),2*(a3+h*z)*(h*a1+x),1+2*h*(1.1+(a1+h*x)*(a2+h*y))];

end