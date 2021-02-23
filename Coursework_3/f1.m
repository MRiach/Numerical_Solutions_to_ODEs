%f_n used in the Lorenz system in exercise 5 and 6 
function fn = f1(x,y,z)

fn = zeros(1,3);
fn(1,1) = 10*(y-x);
fn(1,2) = x*(28-z)-y;
fn(1,3) = x*y-8/3*z;


end