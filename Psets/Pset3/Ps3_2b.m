%Alan Wang
%16.90 Pset #3 problem 2b
%Applying Backward Euler Newton Raphson to Nonlinear Pendulum

%establish time step, iterations, and tolerance
t0 = 0;
tf = 10;
deltat = .01;
n = tf/deltat;
tol = 1e-12;
maxIter = 10;

%set initial conditions and parameters
u0(1) = 0;
u0(2) = pi/4;
g = 9.8;
l = 1;

%establish f(u)
f1 = @(t,u) (-g/l)*sin(u(2));
f2 = @(t,u) u(1);
f = {f1,f2};

%initialize time and solution vectors
t = [t0 zeros(1,n)];
y = [u0.' zeros(2,n)];

%compute Backwards Euler and Newton Raphson
%check with tolerance
for i = 1:n
    error = 1;
    t(i+1) = t(i) + deltat;
    ynew(1) = y(1,i) + deltat*(f1(t(i),y(:,i)));
    ynew(2) = y(2,i) + deltat*(f2(t(i),y(:,i)));
    y(1,i+1) = y(1,i) + deltat*(f1(t(i+1),ynew));
    y(2,i+1) = y(2,i) + deltat*(f2(t(i+1),ynew));
    w = y(:,i);
    res = w - y(:,i) - deltat*[f1(t(i),w);f2(t(i),w)];
    iter = 1;
    while error >= tol
        J = [1, (deltat*g/l)*cos(w(2));-1,1];
        deltaw = J\(-res);
        w = w + deltaw;
        res = w - y(:,i) - deltat*[f1(t(i),w);f2(t(i),w)];
        error = max(abs(res));
        iter = iter + 1;
        if iter == maxIter
            error = 1e-13;
        end
    end
    y(:,i+1) = w;
end

%plot
plot(t,y(2,:))
xlabel('t')
ylabel('theta')
title('Backward Euler Newton Raphson Nonlinear Pendulum')