function dRdw = dRdwfun(Un, deltax)
%compute dRdw for backDiff2

dfdw = @dfdufun;
J = dfdw(Un);
n = size(J);
dRdw = eye(n(1)) - (2/3)*deltax*J;


end

