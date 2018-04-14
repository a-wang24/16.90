function [unew] = NewtonRhapsonBDF2(U,x,dx)

u = U(1,:)';        % u(n)
ulast = U(2,:)';    % u(n-1)

w = u;              % Initial guess: w(n+1) = u(n)
R = -Mfun(w,x)*w + 4/3*Mfun(w,x)*u - 1/3*Mfun(w,x)*ulast + 2/3*dx*Rfun(w,x);

while max(abs(R)) > 1e-6
    dRdw = -dMwdwfun(w,x) + 4/3*dMudwfun(w,u,x) - 1/3*dMudwfun(w,ulast,x) + 2/3*dx*dRdwfun(w,x);
    w = w - (dRdw\R); % dw = -(dRdw)^-1 * R
    R = -Mfun(w,x)*w + 4/3*Mfun(w,x)*u - 1/3*Mfun(w,x)*ulast + 2/3*dx*Rfun(w,x);
end

unew = w; % u(n+1)