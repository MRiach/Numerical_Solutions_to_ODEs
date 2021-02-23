%f_n used in the Robertson chemical reaction in exercise 4
function fn = f4(x,y,z)

fn = zeros(1,3);
fn(1,1) = -0.04*x+10^4*y*z;
fn(1,2) = 0.04*x-10^4*y*z-3*10^7*y^2;
fn(1,3) = 3*10^7*y^2;


end