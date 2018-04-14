%Alan Wang
%16.90 Pset #3 Problem 1b
%Implement Backwards Euler for Temperature Evolution

%set time steps
deltaT = 0.01;
t0 = 0;
tf = 10;
n = tf/deltaT;

%set initial condition and parameter
u0 = 1;
lambda = 0.5;

%set function f
f=@(t,y) -lambda*y;

%initialize time vector and solution vector
t = [t0 zeros(1,n)];
y = [u0 zeros(1,n)];

%compute Backwards Euler numerical solutions
for i = 1:n
    t(i+1) = t(i)+deltaT;
    ynew = y(i) + deltaT*(f(t(i),y(i)));
    y(i+1) = y(i) + deltaT*f(t(i+1),ynew);
end

%find exact solution
t_exact = [t0:deltaT:tf];
y_exact = exp(-lambda*t_exact);

%plot
plot(t_exact,y_exact,t,y)
xlabel('t')
ylabel('u')
title('deltat = 0.01')
legend('Analytical','Backward Euler')