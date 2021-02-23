%f_n used in the RF system in exercise 5
function fn = f5(x,y,z)

fn = zeros(1,3);
fn(1,1) = y*(z-1+x^2)+0.87*x;
fn(1,2) = x*(3*z+1-x^2)+0.87*y;
fn(1,3) = -2*z*(1.1+x*y);


end