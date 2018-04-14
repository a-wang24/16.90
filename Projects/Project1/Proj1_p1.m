%Alan Wang
%16.90 Project 1 Problem 1

%% Setup
% Initialize
clear all
close all

% Parameters
% Declare global variables for supporting functions
global M_hh M_ha M_aa M_ah D_h D_a K_h K_a K_nl
global Q
M_hh = 1; M_ha = .625; M_aa = 1.25; M_ah = 0.25;
D_h = 0.1; D_a = 0.25; %[s^-1]
K_h = 0.2; K_a = 1.25; %[s^-2]
K_nl = 10;

% Initial Conditions
a0 = 0.08; %pitch at t=0
h0 = 0;    %plunge at t=0
p0 = 0;    %da/dt at t=0
v0 = 0;    %dh/dt at t=0
U0 = [0.08 0 0 0]; %initial state vector

% Governing Equations
L = @Lfun;
M = @Mfun;
f1 = @(t,U) U(3);
f2 = @(t,U) U(4);
f3 = @(t,U) (M_ah*D_h*U(4) + M_ah*K_h*U(2) + M_ah*L(U) - M_hh*D_a*U(3) - K_a*M_hh*(1+K_nl*(U(2)^2))*U(1) - M(U)*M_hh)/(M_aa*M_hh - M_ah*M_ha);
f4 = @(t,U) (M_aa*D_h*U(4) + M_aa*K_h*U(2) + M_aa*L(U) - M_ha*D_a*U(3) - K_a*M_ha*(1+K_nl*(U(2)^2))*U(1) - M(U)*M_ha)/(M_ah*M_ha - M_aa*M_hh);
f = {f1,f2,f3,f4};

%Simulation Parameters
t0 = 0;
tf = 60;
deltat = .01;

%% Part 1(i)
% Forward Euler
Q = 1;
[t,y] = forwardEulerMD(t0,tf,U0,deltat,f);
figure(1)
plot(t,y(1:2,:))
legend('a','h')
title('Forward Euler')
xlabel('time [s]')
ylabel('amplitude')

Q = 1.5;
[t,y] = forwardEulerMD(t0,tf,U0,deltat,f);
figure(2)
plot(t,y(1:2,:))
legend('a','h')
title('Forward Euler')
xlabel('time [s]')
ylabel('amplitude')

%% Part 1(ii)
% Midpoint Method
Q = 1;
[t,y] = midpointRuleMD(t0,tf,U0,deltat,f);
figure(3)
plot(t,y(1:2,:))
legend('a','h')
title('Midpoint Method')
xlabel('time [s]')
ylabel('amplitude')

Q = 1.5;
[t,y] = midpointRuleMD(t0,tf,U0,deltat,f);
figure(4)
plot(t,y(1:2,:))
legend('a','h')
title('Midpoint Method')
xlabel('time [s]')
ylabel('amplitude')

%% Part 1(iii)
% Second-Order Backwards Differentiation Scheme (BDF-2)

