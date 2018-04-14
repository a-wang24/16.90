function u2 = u2fun(x)
%u2 piecewise function from 16.90 pset 7
%   Detailed explanation goes here
n = length(x);
u2 = zeros(1,n);
for i = 1:n
    if x(i)>=0 && x(i)<0.5
        u2(i) = exp(x(i));
    elseif x(i)>0.5 && x(i)<=1
        u2(i) = exp(2*(x(i)-.25));
    end
end
end

