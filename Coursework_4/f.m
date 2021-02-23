%f_n used in the IVP in exercise 1 
function fn = f(u,v)

fn = zeros(1,2);
fn(1,1) = -2*u+v;
fn(1,2) = -2*u+v+1e17*(u-v);
end