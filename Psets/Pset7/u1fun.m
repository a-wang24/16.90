function u1 = u1fun(x)
%u1 piecewise function from 16.90 pset 7
%   Detailed explanation goes here

n = length(x);
u1 = zeros(1,n);
for i = 1:n
    if x(i)>=0 && x(i)<0.5
        u1(i) = 5;
    elseif x(i)>0.5 && x(i)<=1
        u1(i) = 1;
    end
end
end

