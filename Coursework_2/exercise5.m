% This function outputs the numbers used in exercise 5
function exercise5()

Y = AB3GlobalError(0,20,[2,3],0.0005);
X = ['The global error for the AB3 method at h = 0.0005 is ',num2str(Y(1))];
disp(X)
Y = AB3GlobalError(0,20,[2,3],0.00025);
X = ['The global error for the AB3 method at h = 0.00025 is ',num2str(Y(1))];
disp(X)
Y = AB3GlobalError(0,20,[2,3],0.000125);
X = ['The global error for the AB3 method at h = 0.000125 is ',num2str(Y(1))];
disp(X)

end 