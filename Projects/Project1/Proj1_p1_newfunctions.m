%Alan Wang
%16.90 Project 1 Problem 1

%% Setup
% Initialize
clear all
close all
warning('off','all');

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
U0 = [a0 h0 p0 v0]; %initial state vector

%Simulation Parameters
t0 = 0;
tf = 60;
deltatf = .01;
deltatm = .0001;

%% Part 1
% initial alpha of 0.08
% numerically obtain solutions for forward Euler, Midpoint and BDF2
% Plot pitch and plunge vs. time for each method
% Start with Q = 1.0 then Q = 1.5
Q = 1.0;
[t1,yf] = forwardEulerMD_functionhandles(t0,tf,U0,deltatf,@fufun);
figure(1)
alpha = subplot(2,1,1);
plunge = subplot(2,1,2);
plot(alpha,t1,yf(1,:))
plot(plunge,t1,yf(2,:))
title(alpha,'Forward Euler Q = 1.0 timestep = .01s')
xlabel('time [s]')
ylabel(alpha, 'alpha')
ylabel(plunge,'h')

[tm,ym] = midpointRuleMD_functionhandles(t0,tf,U0,deltatm,@fufun);
figure(2)
alpha = subplot(2,1,1);
plunge = subplot(2,1,2);
plot(alpha,tm,ym(1,:));
plot(plunge,tm,ym(2,:));
title(alpha,'Midpoint Q = 1.0 timestep = .0001s')
xlabel('time [s]')
ylabel(alpha, 'alpha')
ylabel(plunge,'h')

[tb,yb] = backDiff2(t0,tf,U0,deltatf,@fufun,100,1e-12);
figure(3)
alpha = subplot(2,1,1);
plunge = subplot(2,1,2);
plot(alpha,tb,yb(1,:));
plot(plunge,tb,yb(2,:));
title(alpha,'BDF2 Q = 1.0 timestep = .01s')
xlabel('time [s]')
ylabel(alpha, 'alpha')
ylabel(plunge,'h')

Q=1.5;
[tb,yb] = backDiff2(t0,tf,U0,deltatf,@fufun,100,1e-12);
figure(4)
alpha = subplot(2,1,1);
plunge = subplot(2,1,2);
plot(alpha,tb,yb(1,:));
plot(plunge,tb,yb(2,:));
title(alpha,'BDF2 Q = 1.5 timestep = .01s')
xlabel('time [s]')
ylabel(alpha, 'alpha')
ylabel(plunge,'h')

[t1,y1] = forwardEulerMD_functionhandles(t0,tf,U0,deltatf,@fufun);
figure(5)
alpha = subplot(2,1,1);
plunge = subplot(2,1,2);
plot(alpha,t1,y1(1,:))
plot(plunge,t1,y1(2,:))
title(alpha,'Forward Euler Q = 1.5 timestep = .01s')
xlabel('time [s]')
ylabel(alpha, 'alpha')
ylabel(plunge,'h')

[tm,ym] = midpointRuleMD_functionhandles(t0,tf,U0,deltatm,@fufun);
figure(6)
alpha = subplot(2,1,1);
plunge = subplot(2,1,2);
plot(alpha,tm,ym(1,:));
plot(plunge,tm,ym(2,:));
title(alpha,'Midpoint Q = 1.5 timestep = .0001s')
xlabel('time [s]')
ylabel(alpha, 'alpha')
ylabel(plunge,'h')

%% Part 2
% Find max absolute value of pitch and plunge for all three methods over
% time [0,60]
% vary initial pitch between 0 and 0.08 radians
% Again, compute for Q = 1.0 and Q = 1.5