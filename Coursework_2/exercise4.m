% This function outputs the numbers used in exercise 4 
function exercise4()

Y = ImplicitGlobalError(0,20,[2,3],0.0002);
X = ['The global error for the implicit method at h = 0.0002 is ',num2str(Y(1))];
disp(X)
Y = ImplicitGlobalError(0,20,[2,3],0.0001);
X = ['The global error for the implicit method at h = 0.0001 is ',num2str(Y(1))];
disp(X)
Y = ExplicitGlobalError(0,20,[2,3],0.0001);
X = ['The global error for the explicit method at h = 0.0001 is ',num2str(Y(1))];
disp(X)
Y = ExplicitGlobalError(0,20,[2,3],0.00005);
X = ['The global error for the explicit method at h = 0.00005 is ',num2str(Y(1))];
disp(X)
Y = ExplicitGlobalError(0,20,[2,3],0.00001);
X = ['The global error for the explicit method at h = 0.00001 is ',num2str(Y(1))];
disp(X)

end 