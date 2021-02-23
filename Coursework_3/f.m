%f_n used in the IVP in exercise 2 
function fn = f(u,v,t)

fn = zeros(1,2);
fn(1,1) = -2*u+v+2*sin(t);
fn(1,2) = 998*u-999*v+999*(cos(t)-sin(t));



end