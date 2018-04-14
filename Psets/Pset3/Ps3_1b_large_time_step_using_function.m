%Alan Wang
%16.90 Pset #3 problem 2b
%Applying Backward Euler Newton Raphson to Nonlinear Pendulum using
%written Backward Euler function

clear all
close all

%set parameters
deltaT = 0.01;
t0 = 0;
tf = 10;
y0 = [0 pi/4];
g = 9.8; %m/s^2
l = 1; %m
f = {@(t,y) -(g/l)*sin(u(2)), @(t,y) u(1)};
maxIter = 10;
tol = 1e-12;

%compute Backwards Euler and Newton Raphson
[t,yout] = backwardsEulerMD(t0,tf,y0,deltaT,f);

plot(t,yout(2,:))
xlabel('t')
ylabel('theta')
title('Backward Euler Newton Raphson Nonlinear Pendulum')