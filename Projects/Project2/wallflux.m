% This routine calculates the flux for the Euler 
% equations at a solid wall.
%
% The inputs are:
%
%    UL: state vector in left cell (in the domain)
%     n: normal pointing from boundary into the left cell
% gamma: Ratio of specific heats
%
% The outputs are:
%     H: the flux into the left cell
%  smax: the maximum propagation speed of disturbance
%

function [H, smax] = wallflux(UL, n, gamma)

rL  = UL(1);
uL  = UL(2)/rL;
vL  = UL(3)/rL;
unL = uL*n(1) + vL*n(2);
qL  = sqrt(UL(2)^2 + UL(3)^2)/rL;
utL = sqrt(qL^2 - unL^2);
pL = (gamma-1)*(UL(4) - 0.5*rL*utL^2);
rHL = UL(4) + pL;
cL = sqrt(gamma*pL/rL);

smax = abs(unL) + cL;

H = zeros(size(UL));
H(2) = pL*n(1);
H(3) = pL*n(2);





