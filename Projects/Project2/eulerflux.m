% This routine calculates the flux for the Euler 
% equations using an upwind flux function (specifically
% the HLLE flux function).
%
% The inputs are:
%
%    UL: state vector in left cell
%    UR: state vector in right cell
%     n: normal pointing from the right cell to the left cell
% gamma: Ratio of specific heats
%
% The outputs are:
%     H: the flux into the left cell
%  smax: the maximum propagation speed of disturbance
%

function [H, smax] = eulerflux(UL, UR, n, gamma)

rL = UL(1);
uL = UL(2)/rL;
vL = UL(3)/rL;
unL = uL*n(1) + vL*n(2);
qL = sqrt(UL(2)^2 + UL(3)^2)/rL;
pL = (gamma-1)*(UL(4) - 0.5*rL*qL^2);
rHL = UL(4) + pL;
cL = sqrt(gamma*pL/rL);
sLmin = min(0, unL - cL);
sLmax = max(0, unL + cL);
sLmag = abs(unL) + cL;

HL(1) = rL*unL;
HL(2) = UL(2)*unL + pL*n(1);
HL(3) = UL(3)*unL + pL*n(2);
HL(4) = rHL*unL;

rR = UR(1);
uR = UR(2)/rR;
vR = UR(3)/rR;
unR = uR*n(1) + vR*n(2);
qR = sqrt(UR(2)^2 + UR(3)^2)/rR;
pR = (gamma-1)*(UR(4) - 0.5*rR*qR^2);
rHR = UR(4) + pR;
cR = sqrt(gamma*pR/rR);
sRmin = min(0, unR - cR);
sRmax = max(0, unR + cR);
sRmag = abs(unR) + cR;

HR(1) = rR*unR;
HR(2) = UR(2)*unR + pR*n(1);
HR(3) = UR(3)*unR + pR*n(2);
HR(4) = rHR*unR;

smax = max(sLmag, sRmag);

sLRmin = min( sLmin, sRmin );
sLRmax = max( sLmax, sRmax );

H = 0.5*(HL' + HR') ...
  - 0.5*(sLRmax+sLRmin)/(sLRmax - sLRmin)*(HL' - HR') ...
  + sLRmax*sLRmin/(sLRmax - sLRmin)*(UL-UR);

