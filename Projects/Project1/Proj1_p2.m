%Alan Wang
%16.90 Project 1 Problem 2

%% Part 1
% Use Forward Euler in time and Backward in space to solve

clear all
close all

%import data
Cx = zeros(514);
Cy = zeros(514);
Cx(2:514,2:514) = importdata('Cx.txt');
Cy(2:514,2:514) = importdata('Cy.txt');

%Number of points
Nx = 513;
x = linspace(0,1,Nx+1);
dx = 1/Nx;

Ny = 513;
y = linspace(0,1,Ny+1);
dy = 1/Nx;

% Set initial condition
U0 = zeros(514);

%time step
deltat = 1/1024;
tf = 1;
t0 = 0;
numSteps = (tf-t0)/deltat;

U = cell([1 numSteps+1]);
U{1,1} = U0;

t = [t0 zeros(1,numSteps)];

%Not sure if Cx.txt and Cy.txt are indexed like the images in pset?
%Is Cx(1,1) the bottom left of the image?
%Or is it the top left just like the location of the matrix?
for n = 1:numSteps
    disp(['run: ', num2str(n)]);
    t(n) = n*deltat;
    U{1,n+1} = zeros(514);
    U{1,n+1}(1,:) = sin(t(n)*pi);
    U{1,n+1}(:,1) = sin(t(n)*pi);
    U{1,n+1}(2:514,2:514) = U{1,n}(2:514,2:514) - deltat*(((Cx(2:514,2:514)-Cx(1:513,2:514))/dx).*U{1,n}(2:514,2:514) +...
        Cx(2:514,2:514).*((U{1,n}(2:514,2:514)-U{1,n}(1:513,2:514)))/dx + ...
        ((Cy(2:514,2:514)-Cy(2:514,1:513))/dy).*U{1,n}(2:514,2:514) + ...
        Cy(2:514,2:514).*((U{1,n}(2:514,2:514)-U{1,n}(2:514,1:513))/dy));
%     for i = 1:(Nx+1)
%         for j = 1:(Ny+1)
%             %boundary statement conditional
%             if i == 1 || j == 1
%                 U{1,n}(i,j) = sin(n*pi);
%             else
%                 U{1,n}(i,j) = U{1,n-1}(i,j) - deltat*((((Cx(i,j)-Cx(i-1,j))/dx)*U{1,n-1}(i,j)) ...
%                     + Cx(i,j)*((U{1,n-1}(i,j)-U{1,n-1}(i-1,j))/dx) + ...
%                     ((Cy(i,j)-Cy(i,j-1))/dy)*U{1,n-1}(i,j) + Cy(i,j)*((U{1,n-1}(i,j)-U{1,n-1}(i,j-1))/dy));
%             end
%         end
%     end
end

t1 = .25*1024;     
t2 = .5*1024;
t3 = 1025;

figure(1)
surf(x,x, U{1,t1+1});
shading interp;
caxis([-2 2]);
view(2);
axis equal;
axis([0,1,0,1])
colorbar;
drawnow;
title('t=0.25')
xlabel('x')
ylabel('y')

figure(2)
surf(x,x, U{1,t2+1});
shading interp;
caxis([-2 2]);
view(2);
axis equal;
axis([0,1,0,1])
colorbar;
drawnow;
title('t=0.50')
xlabel('x')
ylabel('y')

figure(3)
surf(x,x, U{1,t3});
shading interp;
caxis([-2 2]);
view(2);
axis equal;
axis([0,1,0,1])
colorbar;
drawnow;
title('t=1')
xlabel('x')
ylabel('y')

