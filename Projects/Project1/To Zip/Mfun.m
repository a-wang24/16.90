function M = Mfun(Un)
%Moment Function
% 
global Q

a = Un(1);

M = -0.7*Q*a;

end

