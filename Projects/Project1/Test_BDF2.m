%test BDF2

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
deltat = .01;

Q = 1;

[y,t] = BDF2(@fufunvect, U0,t0,tf,deltat);
plot(t,y(1:2,:))