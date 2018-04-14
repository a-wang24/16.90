function [x, y] = forwardEulerMD_functionhandlesm(x_i,x_f,y_i,deltax,f)
%forward Euler MD implements the forward Euler method
%this function can handle 1D or multi-dimensional problems
% Input:
% x_i - initial x
% x_f - final x
% y_i - initial y as a horizontal array or scalar
% deltax - step size
% f - the function for y prime as a cell array of functions
% Output:
% y - computed solution
% x - time steps

%calculate number of steps
numSteps = (x_f - x_i)/deltax;

%initialize arrays
x = [x_i zeros(1,numSteps)];
n = size(y_i);
y = [y_i.' zeros(n(2),numSteps)];

%numeric routine
for i = 1:numSteps
    x(i+1) = x(i) + deltax;
    for iter = 1:n(2)
        f_curr = f{iter};
        y(iter, i+1) = y(iter, i) + deltax*f_curr(x(i),y(:,i));
    end
end
end
    
