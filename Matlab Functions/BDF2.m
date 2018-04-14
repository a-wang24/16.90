function [y, x] = BDF2(dydx, y0, x0, xf, dx)
% BDF2 integration function
%
% Input:
% dydx = Function handle so that dy(n)/dx = dydx(y(n))
% y0   = Vector of initial y state
% x0   = Vector of initial x state
% xf   = Vector of final x state
% dx   = Vector of x increments
%
% Output:
% y  = Matrix of y values from BDF2 integration
%      Rows are each x-step of BDF2 integration
%      Columns are components of y (also y0) vector
%      If there are multiple x-ranges (x0, xf, dx), y will be a cell
%      Cells are indices of x-ranges used
% x  = Nx1 matrix of x values from specified vectors of x0, xf, dx
%      Rows are incremental x-steps
%      If there are multiple x-ranges (x0, xf, dx), x will be a cell
%      Cells are indices of x-ranges used

if length(x0) ~= length(xf) || length(x0) ~= length(dx)
    error('Inconsistant set of x: Ensure the lengths of x0, xf, and dx are the same and try again');
end

y = cell(1,length(x0));
x = y;

for xset = 1:length(x0)
    disp(['BDF2 run ' num2str(xset) ' with step size dx = ' num2str(dx(xset))])
    x{xset} = (x0(xset):dx(xset):xf(xset))';
    N = length(x{xset});
    y{xset} = zeros(N,length(y0));
    y{xset}(1,:) = y0(:);
    x{xset} = (x0(xset):dx(xset):xf(xset))';
    
    for n = 1:N-1
        X = x{xset}(n);
        Xnext = x{xset}(n+1);
        if n == 1            
            y{xset}(n+1,:) = y{xset}(n,:) + ...
                dx(xset) * dydx(y{xset}(n,:),X);
        else
            u = [y{xset}(n,:); y{xset}(n-1,:)];
            y{xset}(n+1,:) = NewtonRhapsonBDF2(u,Xnext,dx(xset));
        end
    end
end

if xset == 1
    y = y{1};
    x = x{1};
end

end