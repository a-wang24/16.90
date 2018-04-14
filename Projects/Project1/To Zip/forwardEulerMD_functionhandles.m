function [x, y] = forwardEulerMD_functionhandles(x0,xf,y0,deltax,dfdx)
%forward Euler MD implements the forward Euler method
%this function can handle 1D or multi-dimensional problems
% Input:
% x0 - initial x
% xf - final x
% y0 - initial y as a horizontal array or scalar
% deltax - step size
% dfdx - function handle for df/dx 
% Output:
% y - computed solution
% x - time steps

%calculate number of steps
numSteps = (xf - x0)/deltax;

%initialize arrays
x = [x0 zeros(1,numSteps)];
n = size(y0);
y = [y0.' zeros(n(2),numSteps)];

%numeric routine
for i = 1:numSteps
    x(i+1) = x(i) + deltax;
    y(:,i+1) = y(:,i) + deltax*dfdx(y(:,i));
end
end
    
