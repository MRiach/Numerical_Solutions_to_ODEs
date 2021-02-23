%This is the working used to compute the beta coefficients for the implicit
%method.

format rat

A = [1 2 3; 1/2 2 9/2; 1/6 4/3 9/2];

B = [53/20; 197/60; 85/48];

inv(A)*B