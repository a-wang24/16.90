function wnew = NewtonRaphsonBDF2(vn,vn1,deltax,dfdx,tol,maxCount)
%Applies Newton Raphson to 2nd Order Backwards Differentiation
%   Detailed explanation goes here
w = vn;
Res = @Resfun;
dRdw = @dRdwfun;

Residual = Res(w,vn,vn1,deltax,dfdx);
count = 0;
error = max(abs(Residual));

while error >= tol
    J = dRdw(w,deltax);
    deltaw = J\(-Residual);
    w = w + deltaw;
    Residual = Res(w,vn,vn1,deltax,dfdx);
    error = max(abs(Residual));
    count = count + 1;
    if count >= maxCount
        error = tol*1e-2;
    end
end
wnew = w;
end

