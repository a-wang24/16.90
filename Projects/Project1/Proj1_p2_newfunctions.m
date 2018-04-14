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
h0 = 0;    %plunge at t=0
p0 = 0;    %da/dt at t=0
v0 = 0;    %dh/dt at t=0

%Simulation Parameters
t0 = 0;
tf = 60;
deltatf = .01;
deltatm = .0001;

%% Part 2
% Forward Euler
Q = 1.0;
a0 = [0:.002:.08];
n = length(a0);

% maxaf1 = zeros(1,n);
% maxhf1 = zeros(1,n);
% for i = 1:n
%     U0 = [a0(i) h0 p0 v0];
%     [t1,y1] = forwardEulerMD_functionhandles(t0,tf,U0,deltatf,@fufun);
%     maxaf1(i) = max(abs(y1(1,:)));
%     maxhf1(i) = max(abs(y1(2,:)));
% end
% Q = 1.5;
% maxaf2 = zeros(1,n);
% maxhf2 = zeros(1,n);
% for i = 1:n
%     U0 = [a0(i) h0 p0 v0];
%     [t1,y1] = forwardEulerMD_functionhandles(t0,tf,U0,deltatf,@fufun);
%     maxaf2(i) = max(abs(y1(1,:)));
%     maxhf2(i) = max(abs(y1(2,:)));
% end
% figure(1)
% Q1 = subplot(2,1,1);
% Q15 = subplot(2,1,2);
% plot(Q1,a0,maxaf1,a0,maxhf1)
% plot(Q15,a0,maxaf2,a0,maxhf2)
% legend(Q1, 'abs(a)','abs(h)')
% legend(Q15, 'abs(a)','abs(h)')
% title(Q1, 'Forward Q = 1.0')
% title(Q15, 'Forward Q = 1.5')
% xlabel('alpha at t=0')

%Midpoint Method
% Q=1.0;
% a0 = [0:.005:.08];
% n = length(a0);
% maxam1 = zeros(1,n);
% maxhm1 = zeros(1,n);
% for i = 1:n
%     U0 = [a0(i) h0 p0 v0];
%     [t2,y2] = midpointRuleMD_functionhandles(t0,tf,U0,.001,@fufun);
%     maxam1(i) = max(abs(y2(1,:)));
%     maxhm1(i) = max(abs(y2(2,:)));
% end
% Q = 1.5;
% maxam2 = zeros(1,n);
% maxhm2 = zeros(1,n);
% for i = 1:n
%     U0 = [a0(i) h0 p0 v0];
%     [t2,y2] = midpointRuleMD_functionhandles(t0,tf,U0,.001,@fufun);
%     maxam2(i) = max(abs(y2(1,:)));
%     maxhm2(i) = max(abs(y2(2,:)));
% end
% figure(2)
% Q1 = subplot(2,1,1);
% Q15 = subplot(2,1,2);
% plot(Q1,a0,maxam1,a0,maxhm1)
% plot(Q15,a0,maxam2,a0,maxhm2)
% legend(Q1, 'abs(a)','abs(h)')
% legend(Q15, 'abs(a)','abs(h)')
% title(Q1, 'Midpoint Q = 1.0')
% title(Q15, 'Midpoint Q = 1.5')
% xlabel('alpha at t=0')

%BDF2
maxab1 = zeros(1,n);
maxhb1 = zeros(1,n);
for i = 1:n
    disp(['Q = 1 run: ', num2str(i)]);  
    U0 = [a0(i) h0 p0 v0];
    [t3,y3] = backDiff2(t0,tf,U0,deltatf,@fufun,100,1e-12);
    maxab1(i) = max(abs(y3(1,:)));
    maxhb1(i) = max(abs(y3(2,:)));
end
Q = 1.5;
maxab2 = zeros(1,n);
maxhb2 = zeros(1,n);
for i = 1:n
    disp(['Q = 1.5 run: ', num2str(i)]);  
    U0 = [a0(i) h0 p0 v0];
    [t3,y3] = backDiff2(t0,tf,U0,deltatf,@fufun,100,1e-12);
    maxab2(i) = max(abs(y3(1,:)));
    maxhb2(i) = max(abs(y3(2,:)));
end
figure(3)
Q1 = subplot(2,1,1);
Q15 = subplot(2,1,2);
plot(Q1,a0,maxab1,a0,maxhb1)
plot(Q15,a0,maxab2,a0,maxhb2)
legend(Q1, 'abs(a)','abs(h)')
legend(Q15, 'abs(a)','abs(h)')
title(Q1, 'BDF2 Q = 1.0')
title(Q15, 'BDF2 Q = 1.5')
xlabel('alpha at t=0')




