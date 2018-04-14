function vx = nodalApprox(fx,xi,xf,n)
%interpolate1D - approximates a function as a linear combination of linear
%nodal basis functions
%   fx - function handle for function interpolating
%   phi- function handle for basis functions: phifun
%   xi - initial x
%   xf - final x
%   n  - number of nodal points

xj = linspace(xi,xf,n);
dx = xj(2)-xj(1);
phi = @phifun;
phiCell = phi(xj,dx);
a = zeros(1,n);
points = 49*n+1;
vx = zeros(1,points);

for i = 1:n
    a(i) = fx(xj(i));
    vx = vx + a(i)*phiCell{i};
end
end

