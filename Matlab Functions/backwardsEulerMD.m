function [x, y] = backwardsEulerMD(x_i,x_f,y_i,deltax,f)
% backwards Euler MD implements the backwards Euler method
% this function can handle 1D or multi-dimensional problems
% function f represents y prime
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

%Numeric Routine
%if n is 1, then 1D problem
if n(2) == 1
    f_curr = f{1};
    for i = 1:numSteps
        x(i+1) = x(i)+deltax;
        ynew = y(i)+deltax*(f_curr(x(i),y(i)));
        y(i+1) = y(i)+deltax*f_curr(x(i+1),ynew);
    end
    %if n is not 1, then multi-dimensional problem
else
    for i = 1:numSteps
        x(i+1) = x(i) + deltax;
        ynew = zeros(1,n(2));
        for iter = 1:n(2)
            f_curr = f{iter};
            ynew(iter) = y(iter,i) + deltax*(f_curr(x(i),y(:,i)));
        end
        for iter = 1:n(2)
            f_curr = f{iter};
            y(iter,i+1) = y(iter,i) + deltax*(f_curr(x(i+1),ynew));
        end
    end
end
end

