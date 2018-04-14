function [x, y] = backDiff2(x0,xf,y0,deltax,dfdx,maxCount,tol)
%backDiff2 implements 2nd order Backwards Differentiation 
%using Newton Raphson methods
%this function can handle 1D or multi-dimensional problems
% Input:
% x_i      - initial x
% x_f      - final x
% y_i      - initial y as a horizontal array or scalar
% deltax   - step size
% dfdx     - function handle for df/dx
% maxCount - max iterations to run for Newton Raphson
% Output:
% y - computed solution
% x - time steps

%calculate number of steps
numSteps = (xf - x0)/deltax;

%initialize arrays
x = [x0 zeros(1,numSteps)];
n = size(y0);
y = [y0.' zeros(n(2),numSteps)];

for i = 1:numSteps
    x(i+1) = x(i) + deltax;
    %for first step use forward euler
    if i == 1
        y(:,i+1) = y(:,i) + deltax*dfdx(y(:,i));
    else
        vn = y(:,i);
        vn1 = y(:,i-1);
        y(:,i+1) = NewtonRaphsonBDF2(vn,vn1,deltax,dfdx,tol,maxCount);
    end
end
end

