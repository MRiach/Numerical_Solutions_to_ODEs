%This produces the Jacobian used in Newton method in q4
function matrix = Fprimeq4(x,y,z,a1,a2,a3,h)
 matrix = [1+0.04*h -10^4*h*(a3+h*z) -10^4*h*(a2+h*y);-0.04*h 1+10^4*h*(a3+h*z)+6*10^7*h*(a2+h*y)...
     10^4*h*(a2+h*y);0 -6*10^7*h*(a2+h*y) 1];

end