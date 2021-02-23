%this produces Jacobian used in Newton method
function matrix = Fprime(x,y,z,h,beta)
 matrix = [1+10*h*beta -10*h*beta 0;h*beta*(z-28) 1+h*beta h*beta*x;-h*beta*y -h*beta*x 1+h*beta*8/3];

end
